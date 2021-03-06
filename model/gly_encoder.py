from typing import Optional, Tuple

import torch
import torch.nn as nn
from fairseq.modules import (
    FairseqDropout,
    LayerDropModuleList,
    LayerNorm,
    MultiheadAttention,
    PositionalEmbedding,
    TransformerSentenceEncoderLayer,
)
from fairseq.modules.quant_noise import quant_noise as apply_quant_noise_


def init_bert_params(module):
    """
    Initialize the weights specific to the BERT Model.
    This overrides the default initializations depending on the specified arguments.
        1. If normal_init_linear_weights is set then weights of linear
           layer will be initialized using the normal distribution and
           bais will be set to the specified value.
        2. If normal_init_embed_weights is set then weights of embedding
           layer will be initialized using the normal distribution.
        3. If normal_init_proj_weights is set then weights of
           in_project_weight for MultiHeadAttention initialized using
           the normal distribution (to be validated).
    """

    if isinstance(module, nn.Linear):
        module.weight.data.normal_(mean=0.0, std=0.02)
        if module.bias is not None:
            module.bias.data.zero_()
    if isinstance(module, nn.Embedding):
        module.weight.data.normal_(mean=0.0, std=0.02)
        if module.padding_idx is not None:
            module.weight.data[module.padding_idx].zero_()
    if isinstance(module, MultiheadAttention):
        module.q_proj.weight.data.normal_(mean=0.0, std=0.02)
        module.k_proj.weight.data.normal_(mean=0.0, std=0.02)
        module.v_proj.weight.data.normal_(mean=0.0, std=0.02)


class GlyEncoder(nn.Module):
    """
    Implementation for a Bi-directional Transformer based Sentence Encoder used
    in BERT/XLM style pre-trained models.
    This first computes the token embedding using the token embedding matrix,
    position embeddings (if specified) and segment embeddings
    (if specified). After applying the specified number of
    TransformerEncoderLayers, it outputs all the internal states of the
    encoder as well as the final representation associated with the first
    token (usually CLS token).
    Input:
        - tokens: B x T matrix representing sentences
        - segment_labels: B x T matrix representing segment label for tokens
    Output:
        - a tuple of the following:
            - a list of internal model states used to compute the
              predictions where each tensor has shape T x B x C
            - sentence representation associated with first input token
              in format B x C.
    """

    def __init__(
        self,
        padding_idx: int,
        vocab_size: int,
        num_encoder_layers: int = 6,
        embedding_dim: int = 768,
        ffn_embedding_dim: int = 3072,
        num_attention_heads: int = 8,
        dropout: float = 0.1,
        attention_dropout: float = 0.1,
        activation_dropout: float = 0.1,
        layerdrop: float = 0.0,
        max_seq_len: int = 256,
        num_branch: int = 2**10,
        max_branch: int = 50,
        connections: int = 2**9,
        aromatic: int = 3,
        parent=10,
        encoder_normalize_before: bool = False,
        apply_bert_init: bool = False,
        activation_fn: str = "relu",
        learned_pos_embedding: bool = True,
        embed_scale: float = None,
        freeze_embeddings: bool = False,
        n_trans_layers_to_freeze: int = 0,
        export: bool = False,
        traceable: bool = False,
        q_noise: float = 0.0,
        qn_block_size: int = 8,
    ) -> None:

        super().__init__()
        self.gen=False
        
        self.padding_idx = padding_idx
        self.vocab_size = vocab_size
        self.dropout_module = FairseqDropout(
            dropout, module_name=self.__class__.__name__
        )
        self.layerdrop = layerdrop
        self.max_seq_len = max_seq_len
        self.embedding_dim = embedding_dim
        self.num_branch = num_branch
        self.max_branch = max_branch
        self.connections = connections
        self.aromatic = aromatic
        self.parent=parent
        self.apply_bert_init = apply_bert_init
        self.traceable = traceable

        self.embed_tokens = nn.Embedding(self.vocab_size, self.embedding_dim, self.padding_idx)
        
        self.embed_scale = embed_scale

        if q_noise > 0:
            self.quant_noise = apply_quant_noise_(
                nn.Linear(self.embedding_dim, self.embedding_dim, bias=False),
                q_noise,
                qn_block_size,
            )
        else:
            self.quant_noise = None

        self.aromatic_embeddings = nn.Embedding(self.aromatic, self.embedding_dim, padding_idx=self.padding_idx)
        self.branch_embeddings = nn.Embedding(self.num_branch, self.embedding_dim, padding_idx=self.padding_idx)
        self.depth_embeddings = nn.Embedding(self.max_branch, self.embedding_dim)

        self.connection_embeddings = nn.Embedding(self.connections, self.embedding_dim, padding_idx=self.padding_idx)
        self.parent_embeddings = nn.Embedding(self.parent, self.embedding_dim, padding_idx=self.padding_idx)
        if self.layerdrop > 0.0:
            self.layers = LayerDropModuleList(p=self.layerdrop)
        else:
            self.layers = nn.ModuleList([])
        self.layers.extend(
            [
                self.build_transformer_sentence_encoder_layer(
                    embedding_dim=self.embedding_dim,
                    ffn_embedding_dim=ffn_embedding_dim,
                    num_attention_heads=num_attention_heads,
                    dropout=self.dropout_module.p,
                    attention_dropout=attention_dropout,
                    activation_dropout=activation_dropout,
                    activation_fn=activation_fn,
                    export=export,
                    q_noise=q_noise,
                    qn_block_size=qn_block_size,
                )
                for _ in range(num_encoder_layers)
            ]
        )

        if encoder_normalize_before:
            self.emb_layer_norm = LayerNorm(self.embedding_dim, export=export)
        else:
            self.emb_layer_norm = None

        # Apply initialization of model params after building the model
        if self.apply_bert_init:
            self.apply(init_bert_params)

        def freeze_module_params(m):
            if m is not None:
                for p in m.parameters():
                    p.requires_grad = False

        if freeze_embeddings:
            freeze_module_params(self.embed_tokens)
            freeze_module_params(self.segment_embeddings)
            freeze_module_params(self.embed_positions)
            freeze_module_params(self.emb_layer_norm)

        for layer in range(n_trans_layers_to_freeze):
            freeze_module_params(self.layers[layer])


    def build_transformer_sentence_encoder_layer(
        self,
        embedding_dim,
        ffn_embedding_dim,
        num_attention_heads,
        dropout,
        attention_dropout,
        activation_dropout,
        activation_fn,
        export,
        q_noise,
        qn_block_size,
    ):
        return TransformerSentenceEncoderLayer(
            embedding_dim=embedding_dim,
            ffn_embedding_dim=ffn_embedding_dim,
            num_attention_heads=num_attention_heads,
            dropout=dropout,
            attention_dropout=attention_dropout,
            activation_dropout=activation_dropout,
            activation_fn=activation_fn,
            export=export,
            q_noise=q_noise,
            qn_block_size=qn_block_size,
        )

    def forward(
        self,
        allx: torch.Tensor,
        last_state_only: bool = True,
        token_embeddings: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        is_tpu = allx.device.type == "xla"

        # compute padding mask. This is needed for multi-head attention
        
        
        if not self.gen:
            padding_mask = allx[:,:,0].eq(self.padding_idx)
            if not self.traceable and not is_tpu and not padding_mask.any():
                padding_mask = None
        else:
            padding_mask = None
            
        

        if token_embeddings is not None:
            x = token_embeddings
        else:
            if not self.gen:
                x = self.embed_tokens(allx[:,:,0])
            else:
                x= self.embed_tokens(allx)
#                 print(x.size())
                
#         print(x.size())
        
        if self.embed_scale is not None:
            x = x * self.embed_scale

        if self.connection_embeddings is not None:
            x = x + self.connection_embeddings(allx[:,:,2])

        if self.branch_embeddings is not None:
            x = x + self.branch_embeddings(allx[:,:,3])
            
        if self.depth_embeddings is not None:
            x = x + self.depth_embeddings(allx[:,:,5])
            
        if self.aromatic_embeddings is not None:
            x = x + self.aromatic_embeddings(allx[:,:,1])
            
        if self.parent_embeddings is not None:
            x = x + self.parent_embeddings(allx[:,:,4])

        if self.quant_noise is not None:
            x = self.quant_noise(x)

        if self.emb_layer_norm is not None:
            x = self.emb_layer_norm(x)

#         if not self.gen:
        x = self.dropout_module(x)

        # account for padding while computing the representation
        if padding_mask is not None:
            x = x * (1 - padding_mask.unsqueeze(-1).type_as(x))

        # B x T x C -> T x B x C
        x = x.transpose(0, 1)
        inner_states = []
        if not last_state_only:
            inner_states.append(x)
            
        for layer in self.layers:
            x, _ = layer(x, self_attn_padding_mask=padding_mask)
            if not last_state_only:
                inner_states.append(x)

        sentence_rep = x[0, :, :]

        if last_state_only:
            inner_states = [x]

#         if self.traceable:
#             return torch.stack(inner_states), sentence_rep
#         else:
        return inner_states, sentence_rep


class GlyEncoder_el(nn.Module):

    def __init__(
        self,
        padding_idx: int,
        vocab_size: int,
        num_encoder_layers: int = 6,
        embedding_dim: int = 768,
        ffn_embedding_dim: int = 3072,
        num_attention_heads: int = 8,
        dropout: float = 0.1,
        attention_dropout: float = 0.1,
        activation_dropout: float = 0.1,
        layerdrop: float = 0.0,
        max_seq_len: int = 256,
        num_branch: int = 11,
        max_dep: int = 50,
        connections: int = 9,
        aromatic: int = 3,
        parent: int =10,
        encoder_normalize_before: bool = False,
        apply_bert_init: bool = False,
        activation_fn: str = "relu",
        learned_pos_embedding: bool = True,
        embed_scale: float = None,
        freeze_embeddings: bool = False,
        n_trans_layers_to_freeze: int = 0,
        export: bool = False,
        traceable: bool = False,
        q_noise: float = 0.0,
        qn_block_size: int = 8,
    ) -> None:

        super().__init__()
        self.padding_idx = padding_idx
        self.vocab_size = vocab_size
        self.dropout_module = FairseqDropout(
            dropout, module_name=self.__class__.__name__
        )
        self.layerdrop = layerdrop
        self.max_seq_len = max_seq_len
        self.embedding_dim = embedding_dim
        self.num_branch = num_branch
        self.max_dep = max_dep
        self.connections = connections
        self.aromatic = aromatic
        self.apply_bert_init = apply_bert_init
        self.traceable = traceable
        self.parent=parent
        self.embed_scale = embed_scale

        if q_noise > 0:
            self.quant_noise = apply_quant_noise_(
                nn.Linear(self.embedding_dim, self.embedding_dim, bias=False),
                q_noise,
                qn_block_size,
            )
        else:
            self.quant_noise = None
            
        
        self.embed_tokens = nn.Embedding(self.vocab_size, int(self.embedding_dim/2), self.padding_idx)
        self.aromatic_embeddings = nn.Embedding(self.aromatic, int(self.embedding_dim/4), padding_idx=self.padding_idx)
        self.branch_embeddings = nn.Conv1d(self.num_branch, int(self.embedding_dim/4),1)
        self.depth_embeddings = nn.Embedding(self.max_dep, int(self.embedding_dim/4))
        self.connection_embeddings = nn.Conv1d(self.connections, int(self.embedding_dim/4),1)
        self.parent_embeddings = nn.Embedding(self.parent, int(self.embedding_dim/4), padding_idx=self.padding_idx)
        
        if self.layerdrop > 0.0:
            self.layers = LayerDropModuleList(p=self.layerdrop)
        else:
            self.layers = nn.ModuleList([])
        self.layers.extend(
            [
                self.build_transformer_sentence_encoder_layer(
                    embedding_dim=self.embedding_dim,
                    ffn_embedding_dim=ffn_embedding_dim,
                    num_attention_heads=num_attention_heads,
                    dropout=self.dropout_module.p,
                    attention_dropout=attention_dropout,
                    activation_dropout=activation_dropout,
                    activation_fn=activation_fn,
                    export=export,
                    q_noise=q_noise,
                    qn_block_size=qn_block_size,
                )
                for _ in range(num_encoder_layers)
            ]
        )

        if encoder_normalize_before:
            self.emb_layer_norm = LayerNorm(self.embedding_dim, export=export)
        else:
            self.emb_layer_norm = None

        # Apply initialization of model params after building the model
        if self.apply_bert_init:
            self.apply(init_bert_params)

        def freeze_module_params(m):
            if m is not None:
                for p in m.parameters():
                    p.requires_grad = False

        if freeze_embeddings:
            freeze_module_params(self.embed_tokens)
            freeze_module_params(self.segment_embeddings)
            freeze_module_params(self.embed_positions)
            freeze_module_params(self.emb_layer_norm)

        for layer in range(n_trans_layers_to_freeze):
            freeze_module_params(self.layers[layer])


    def build_transformer_sentence_encoder_layer(
        self,
        embedding_dim,
        ffn_embedding_dim,
        num_attention_heads,
        dropout,
        attention_dropout,
        activation_dropout,
        activation_fn,
        export,
        q_noise,
        qn_block_size,
    ):
        return TransformerSentenceEncoderLayer(
            embedding_dim=embedding_dim,
            ffn_embedding_dim=ffn_embedding_dim,
            num_attention_heads=num_attention_heads,
            dropout=dropout,
            attention_dropout=attention_dropout,
            activation_dropout=activation_dropout,
            activation_fn=activation_fn,
            export=export,
            q_noise=q_noise,
            qn_block_size=qn_block_size,
        )

    def forward(
        self,
        allx: torch.Tensor,
        last_state_only: bool = True,
        token_embeddings: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        is_tpu = allx.device.type == "xla"

        # compute padding mask. This is needed for multi-head attention
        padding_mask = allx[:,:,0].eq(self.padding_idx)
        if not self.traceable and not is_tpu and not padding_mask.any():
            padding_mask = None

        if token_embeddings is not None:
            tx = token_embeddings
        else:
            tx = self.embed_tokens(allx[:,:,0])

        if self.embed_scale is not None:
            tx = tx * self.embed_scale

        if self.connection_embeddings is not None:
            cx = torch.transpose(self.connection_embeddings(torch.transpose(allx[:,:,4:13].float(),1,2)),1,2)

        if self.branch_embeddings is not None:
            bx = torch.transpose(self.branch_embeddings(torch.transpose(allx[:,:,13:].float(),1,2)),1,2)
            
        if self.depth_embeddings is not None:
            bx = bx + self.depth_embeddings(allx[:,:,3])
            
        if self.aromatic_embeddings is not None:
            cx = cx + self.aromatic_embeddings(allx[:,:,1])
            
        if self.parent_embeddings is not None:
            cx = cx + self.parent_embeddings(allx[:,:,2])
            
        x=torch.cat([tx,bx,cx],2)

        if self.quant_noise is not None:
            x = self.quant_noise(x)

        if self.emb_layer_norm is not None:
            x = self.emb_layer_norm(x)

        x = self.dropout_module(x)

        # account for padding while computing the representation
        if padding_mask is not None:
            x = x * (1 - padding_mask.unsqueeze(-1).type_as(x))

        # B x T x C -> T x B x C
        x = x.transpose(0, 1)

        inner_states = []
        if not last_state_only:
            inner_states.append(x)

        for layer in self.layers:
            x, _ = layer(x, self_attn_padding_mask=padding_mask)
            if not last_state_only:
                inner_states.append(x)

        sentence_rep = x[0, :, :]

        if last_state_only:
            inner_states = [x]

#         if self.traceable:
#             return torch.stack(inner_states), sentence_rep
#         else:
        return inner_states, sentence_rep