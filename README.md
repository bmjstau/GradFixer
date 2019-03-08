CTF usually does the coordinate system transform from
dewar -> subject head space when the experimenter runs 
the headlocalization at the beginning of a run.
If the subject's position drifts over the run, this will
lead to source-localised positions being wrong by the amount
of drift, as the MRI-based sourcespace is also aligned to 
the subject head space.

To minimize this error, it would be better to align the
subject head space coordinate system not to the subject's 
initial, but to it's average position during an experiment.
This small toolbox achieves this by transforming the 
gradiometer coordinates of a dataset from the initial to
the average head position. This should be done before
computing the subject-specific source models as well as the
leadfield.