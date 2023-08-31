#!/bin/bash

../../vc_mesh -vc vc.in
python3 check_points_in_mesh.py --num_pts 5
