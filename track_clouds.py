#!/usr/bin/env python
import ConfigParser, sys
import cloudtracker.main

if len(sys.argv) < 2:
    raise "No config file given"

config = ConfigParser.RawConfigParser()

config.read(sys.argv[1])
model_config = {}
for option in ('ug', 'vg', 'dt', 'dz', 'dy', 'dx'):
    model_config[option] = config.getfloat('modelconfig', option)
for option in ('nt', 'nz', 'ny', 'nx'):
    model_config[option] = config.getint('modelconfig', option)
model_config['input_directory'] = config.get('modelconfig', 'input_directory')

cloudtracker.main.main(model_config)
