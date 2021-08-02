#!/bin/bash

~/zarusoba_viewer/tools/calculate_nodal_fields_from_integration_point_fields.py mesh.msh Strain\@IntegrationPoint.dat Strain.dat
~/zarusoba_viewer/tools/calculate_nodal_fields_from_integration_point_fields.py mesh.msh Stress\@IntegrationPoint.dat Stress.dat
~/zarusoba_viewer/src/zarusoba_viewer.exe mesh.msh Stress.dat Strain.dat displacement.dat

