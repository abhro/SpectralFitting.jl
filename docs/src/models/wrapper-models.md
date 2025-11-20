# Wrapper models

Wrapper models are a type of "meta-model" that can be used to change the behaviour of other models. Their constructors commonly take a model as an argument, and change the way this model is called.

The following wrapper models are available

```@index
Pages = ["wrapper-models.md"]
Order = [:type]
```

```@docs
ParameterPatch
apply_patch!
AsConvolution
AutoCache
AbstractModelWrapper
```