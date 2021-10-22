2021/10/22 (Jan)
================

* add `malevnc_v0` source data and ingest parameters
* splits:
  * `skeleton`
    * 0.8 for (train ∪ validation) / test
    * 0.8 for train / validation
    * meant for training, hyperparameter search
  * `skeleton_no_test`
    * 1.0 for (train ∪ validation) / test
    * 0.8 for train / validation
    * meant for production on best `skeleton` network
