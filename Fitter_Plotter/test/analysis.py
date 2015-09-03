model = build_model_from_rootfile('result.root', include_mc_uncertainties=True)
model.fill_histogram_zerobins()

for p in model.distribution.get_parameters():
    model.distribution.set_distribution_parameters(p, range = [-5.0, 5.0])
signal_process_groups = {'': []}
parVals = mle(model, input = 'data', n=1, signal_process_groups = signal_process_groups)
print parVals

parameter_values = {}
for p in model.get_parameters([]):
  parameter_values[p] = parVals[''][p][0][0]

hist = evaluate_prediction(model, parameter_values, include_signal = False)
write_histograms_to_rootfile(hist, 'ttbar_theta_histos.root')

model_summary(model)
report.write_html('htmlout')
