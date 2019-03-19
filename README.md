# LGADS

LGADs is a Python script for analyzing electrical tests of diodes (incl. LGADs) for the HGTD project. It can overlay plots of IV, CV, and doping profile curves for sets of sensors.

## Usage

Developed for python3.5, requires additional packages.

### 1) Specify directory containting all test results across all sites
```
directory = main_directory/
```

### 2) Initialize a site (i.e. collection of sensors across one or more wafers specific to a production site)
* the local folder structure must be /directory/path_to_site/
* provide a string as a name for the site (**NOTE:** an output directory will be created with this name)
* provide a list of strings of all wafers specific to the site
    * the string names provided must correspond to the names of the subdirectories within path_to_site/
    * within each wafer subdirectory there must be separate directories IV/ and CV/ containing all sensor tests specific to the wafer
* provide the area of the sensor in m**2
```
ExampleSite = Site(path_to_site, site_name, list_of_wafers, area_of_sensor)
```

### 3) Plot all IV and CV results
```
ExampleSite.show()
```