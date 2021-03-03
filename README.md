# MUSEROS Version 2020

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

MUSER is a Solor-Dedicated Radiograph that was built in Inner-Mongoria of China. 

The MUSEROS is the data processing application for the MUSER. 

This repository contains:

1. Source codes of the data processing pipeline. 
2. Applications of the MUSEROS.
3. Other examples. 


## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)
	- [Generator](#generator)
- [Badge](#badge)
- [Example Readmes](#example-readmes)
- [Related Efforts](#related-efforts)
- [Maintainers](#maintainers)
- [Contributing](#contributing)
- [License](#license)

## Background

This respository include the 2nd version of the MUSEROS, which supports data processing for the MUSER. It's dependant on the [RASCIL](https://gitlab.com/skatelescope/external/rascil).

## Download

This project uses  [RASCIL](https://gitlab.com/skatelescope/external/rascil) that should be installed in advance. Please verify the environmental variables. Meahwhile, casacore and it's dependancies should be installed. 

```sh
$ git clone https://github.com/astronomical-data-processing/museros2020
```

## Installation

In general, the environmental variables should be setup for the MUSEROS2020. 

```sh
# MUSER - source location 
export MUSER=/home/muser/museros2020
# DATA location
export MUSER_DATA=/opt/archive/MUSER-1/dat
# PYTHONPATH
export PYTHONPATH=/Users/meiying/anaconda3/lib/python3.7/site-packages/:$RASCIL:$MUSER:$PYTHONPATH

```

### Contributors


## License

[MIT](LICENSE) © ASTROLAB