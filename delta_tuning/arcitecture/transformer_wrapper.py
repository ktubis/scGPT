import torch
from torch import nn
from typing import Optional, Tuple
import torch.nn.functional as F
from torch.nn.init import constant_, xavier_uniform_


# This function is used to copy the weights of the model trained with the original
# MultiheadAttention module to the model trained with the deconstructed version of the
# MultiheadAttention module.
def copy_original_model(pretrained_model_dict, model_dict):
    if "transformer_encoder.layers.0.self_attn.Q.weight" in model_dict and "transformer_encoder.layers.0.self_attn.Q.weight" not in pretrained_model_dict:
        i = 0
        while f"transformer_encoder.layers.{i}.self_attn.Q.weight" in model_dict:
            w_q, w_k, w_v = pretrained_model_dict[f"transformer_encoder.layers.{i}.self_attn.in_proj_weight"].chunk(3, dim=0)
            b_q, b_k, b_v = pretrained_model_dict[f"transformer_encoder.layers.{i}.self_attn.in_proj_bias"].chunk(3, dim=0)
            model_dict[f"transformer_encoder.layers.{i}.self_attn.Q.weight"] = w_q
            model_dict[f"transformer_encoder.layers.{i}.self_attn.K.weight"] = w_k
            model_dict[f"transformer_encoder.layers.{i}.self_attn.V.weight"] = w_v
            model_dict[f"transformer_encoder.layers.{i}.self_attn.Q.bias"] = b_q
            model_dict[f"transformer_encoder.layers.{i}.self_attn.K.bias"] = b_k
            model_dict[f"transformer_encoder.layers.{i}.self_attn.V.bias"] = b_v
            i += 1

class MultiheadAttentionDeconstructed(nn.MultiheadAttention):
    """
    This is a deconstructed version of the MultiheadAttention module.
    It explicitely defines the Q, K, V of the MultiheadAttention module.
    It is used for the LoRA delta method that requires the Q, K, V to be of type nn.Module.
    Only the initialization differs, all the other methods are the same as the original MultiheadAttention module.
    """

    __constants__ = ['batch_first']
    bias_k: Optional[torch.Tensor]
    bias_v: Optional[torch.Tensor]

    def __init__(self, embed_dim, num_heads, dropout=0., bias=True, add_bias_kv=False, add_zero_attn=False,
                 kdim=None, vdim=None, batch_first=False, device=None, dtype=None) -> None:
        factory_kwargs = {'device': device, 'dtype': dtype}
        super().__init__(embed_dim, num_heads, dropout, bias, add_bias_kv,
                         add_zero_attn, kdim, vdim, batch_first, **factory_kwargs)

        assert self.kdim == embed_dim and self.vdim == embed_dim
        assert self._qkv_same_embed_dim == True
        assert add_bias_kv is False, "add_bias_kv is not supported in this implementation"

        self.Q = nn.Linear(embed_dim, embed_dim, bias=True, **factory_kwargs)
        self.K = nn.Linear(embed_dim, embed_dim, bias=True, **factory_kwargs)
        self.V = nn.Linear(embed_dim, embed_dim, bias=True, **factory_kwargs)

        xavier_uniform_(self.Q.weight)
        xavier_uniform_(self.K.weight)
        xavier_uniform_(self.V.weight)
        constant_(self.Q.bias, 0.)
        constant_(self.V.bias, 0.)
        constant_(self.K.bias, 0.)
        

    def forward(
            self,
            query: torch.Tensor,
            key: torch.Tensor,
            value: torch.Tensor,
            key_padding_mask: Optional[torch.Tensor] = None,
            need_weights: bool = True,
            attn_mask: Optional[torch.Tensor] = None,
            average_attn_weights: bool = True,
            is_causal : bool = False) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
        r"""
        This function replicates the forward function of the original MultiheadAttention module,
        but uses the Q, K, V defined in the __init__ function as the basis for in_proj_weight and in_proj_bias,
        which are used instead of the original self.in_proj_weight and self.in_proj_bias.
        """

        in_proj_weight = torch.cat([self.Q.weight, self.K.weight, self.V.weight], dim=0)
        in_proj_bias = torch.cat([self.Q.bias, self.K.bias, self.V.bias], dim=0)

        why_not_fast_path = ''
        if ((attn_mask is not None and torch.is_floating_point(attn_mask))
           or (key_padding_mask is not None) and torch.is_floating_point(key_padding_mask)):
            why_not_fast_path = "floating-point masks are not supported for fast path."

        is_batched = query.dim() == 3

        key_padding_mask = F._canonical_mask(
            mask=key_padding_mask,
            mask_name="key_padding_mask",
            other_type=F._none_or_dtype(attn_mask),
            other_name="attn_mask",
            target_type=query.dtype
        )

        attn_mask = F._canonical_mask(
            mask=attn_mask,
            mask_name="attn_mask",
            other_type=None,
            other_name="",
            target_type=query.dtype,
            check_other=False,
        )


        if not is_batched:
            why_not_fast_path = f"input not batched; expected query.dim() of 3 but got {query.dim()}"
        elif query is not key or key is not value:
            # When lifting this restriction, don't forget to either
            # enforce that the dtypes all match or test cases where
            # they don't!
            why_not_fast_path = "non-self attention was used (query, key, and value are not the same Tensor)"
        elif self.in_proj_bias is not None and query.dtype != self.in_proj_bias.dtype:
            why_not_fast_path = f"dtypes of query ({query.dtype}) and self.in_proj_bias ({self.in_proj_bias.dtype}) don't match"
        elif in_proj_weight is None:
            why_not_fast_path = "in_proj_weight was None"
        elif query.dtype != in_proj_weight.dtype:
            # this case will fail anyway, but at least they'll get a useful error message.
            why_not_fast_path = f"dtypes of query ({query.dtype}) and in_proj_weight ({in_proj_weight.dtype}) don't match"
        elif self.training:
            why_not_fast_path = "training is enabled"
        elif (self.num_heads % 2) != 0:
            why_not_fast_path = "self.num_heads is not even"
        elif not self.batch_first:
            why_not_fast_path = "batch_first was not True"
        elif self.bias_k is not None:
            why_not_fast_path = "self.bias_k was not None"
        elif self.bias_v is not None:
            why_not_fast_path = "self.bias_v was not None"
        elif self.add_zero_attn:
            why_not_fast_path = "add_zero_attn was enabled"
        elif not self._qkv_same_embed_dim:
            why_not_fast_path = "_qkv_same_embed_dim was not True"
        elif query.is_nested and (key_padding_mask is not None or attn_mask is not None):
            why_not_fast_path = "supplying both src_key_padding_mask and src_mask at the same time \
                                 is not supported with NestedTensor input"
        elif torch.is_autocast_enabled():
            why_not_fast_path = "autocast is enabled"

        if not why_not_fast_path:
            tensor_args = (
                query,
                key,
                value,
                in_proj_weight,
                in_proj_bias,
                self.out_proj.weight,
                self.out_proj.bias,
            )
            # We have to use list comprehensions below because TorchScript does not support
            # generator expressions.
            if torch.overrides.has_torch_function(tensor_args):
                why_not_fast_path = "some Tensor argument has_torch_function"
            elif nn.modules._is_make_fx_tracing():
                why_not_fast_path = "we are running make_fx tracing"
            elif not all(nn.modules._check_arg_device(x) for x in tensor_args):
                why_not_fast_path = ("some Tensor argument's device is neither one of "
                                     f"cpu, cuda or {torch.utils.backend_registration._privateuse1_backend_name}")
            elif torch.is_grad_enabled() and any(nn.modules._arg_requires_grad(x) for x in tensor_args):
                why_not_fast_path = ("grad is enabled and at least one of query or the "
                                     "input/output projection weights or biases requires_grad")
            if not why_not_fast_path:
                merged_mask, mask_type = self.merge_masks(attn_mask, key_padding_mask, query)

                if in_proj_bias is not None and in_proj_weight is not None:
                    return torch._native_multi_head_attention(
                        query,
                        key,
                        value,
                        self.embed_dim,
                        self.num_heads,
                        in_proj_weight,
                        in_proj_bias,
                        self.out_proj.weight,
                        self.out_proj.bias,
                        merged_mask,
                        need_weights,
                        average_attn_weights,
                        mask_type)

        any_nested = query.is_nested or key.is_nested or value.is_nested
        assert not any_nested, ("MultiheadAttention does not support NestedTensor outside of its fast path. " +
                                f"The fast path was not hit because {why_not_fast_path}")

        if self.batch_first and is_batched:
            # make sure that the transpose op does not affect the "is" property
            if key is value:
                if query is key:
                    query = key = value = query.transpose(1, 0)
                else:
                    query, key = (x.transpose(1, 0) for x in (query, key))
                    value = key
            else:
                query, key, value = (x.transpose(1, 0) for x in (query, key, value))


        attn_output, attn_output_weights = F.multi_head_attention_forward(
            query, key, value, self.embed_dim, self.num_heads,
            in_proj_weight, in_proj_bias,
            self.bias_k, self.bias_v, self.add_zero_attn,
            self.dropout, self.out_proj.weight, self.out_proj.bias,
            training=self.training,
            key_padding_mask=key_padding_mask,
            need_weights=need_weights,
            attn_mask=attn_mask,
            average_attn_weights=average_attn_weights,
            is_causal=is_causal)
        if self.batch_first and is_batched:
            return attn_output.transpose(1, 0), attn_output_weights
        else:
            return attn_output, attn_output_weights

    
class CustomTransformerEncoderLayer(nn.TransformerEncoderLayer):
    def __init__(self, *args, attention_cls=None, **kwargs):
        super().__init__(*args, **kwargs)
        if attention_cls is not None:
            self.self_attn = attention_cls(
                self.self_attn.embed_dim,
                self.self_attn.num_heads,
                dropout=self.self_attn.dropout,
                bias=True,
                batch_first=self.self_attn.batch_first,
                device=self.self_attn.in_proj_weight.device,
                dtype=self.self_attn.in_proj_weight.dtype,
            )