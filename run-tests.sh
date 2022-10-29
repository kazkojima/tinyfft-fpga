#!/bin/bash

export GENERATE_VCDS=1

python3 -m unittest vfft.vfft.FixedPointFFTTest
