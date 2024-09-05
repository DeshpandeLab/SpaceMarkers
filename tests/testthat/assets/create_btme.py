from bayestme import data, synthetic_data
stdata = synthetic_data.generate_demo_dataset()

from bayestme import deconvolution
from bayestme.common import InferenceType
best_spatial_smoothing_parameter = 1000.0
best_n_components = 3

deconvolution_result = deconvolution.sample_from_posterior(
    data=stdata,
    n_components=best_n_components,
    spatial_smoothing_parameter=best_spatial_smoothing_parameter,
    n_samples=100,
    n_svi_steps=10_000,
    expression_truth=None,
    use_spatial_guide=True)

data.add_deconvolution_results_to_dataset(
    stdata, deconvolution_result
)

stdata.save('btme.h5ad')