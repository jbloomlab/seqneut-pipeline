# CHANGELOG

## version 1.0.0
- Separate this pipeline from the configuration and code for specific analyses. This is important because we want to be able to update the pipeline in a way that can be used across multiple analyses. So the configuration and results are no longer part of the pipeline, so they can be configured for individual studies while using the same base code.
- Rename from names like "NGS neuts" to "seqneut", because in paper we are choosing to describe as sequencing-based neutralization assays.
- Use `neutcurve` 0.6.0, which is updated to be easier to maintain 
