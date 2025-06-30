# This folder should contain sub-folders with user specific configurations.

For instance, to create a specific user configuration named 'my_config':

* create a folder named 'my_config' in 'configs' folder

* place all your user specific models and stock files inside 'my_config' folder

* place config.yml file with configurations next to models inside 'my_config' file. Paths to models should be recorded as follows:
 /app/configs/<your_custom_config_name>/uspto_model.onnx
Proper paths naming is very important!!! Below is the example of config.yml contents:

search:
  max_transforms: 3
expansion:
  uspto:
    - /app/configs/my_config/uspto_model.onnx
    - /app/aconfigs/my_config/uspto_templates.csv.gz
  ringbreaker:
    - /app/configs/my_config/uspto_ringbreaker_model.onnx
    - /app/configs/my_config/uspto_ringbreaker_templates.csv.gz
filter:
  uspto: /app/configs/my_config/uspto_filter_model.onnx
stock:
  zinc: /app/configs/my_config/zinc_stock.hdf5

Once created, you can select your custom configuration using settings icon on Retrosynthesis widget.
