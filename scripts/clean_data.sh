#!/usr/bin/env bash

[ -f 'TRACK_CHEMISTRY' ]      && rm  TRACK_CHEMISTRY 
[ -f 'OUTPUT' ]               && rm  OUTPUT  
[ -f 'TAGGED_TRAJECTORY' ]    && rm  TAGGED_TRAJECTORY
[ -f 'OCF' ]                  && rm  OCF 
[ -f 'OCF_AVG' ]              && rm  OCF_AVG 
[ -f 'MSD' ]                  && rm  MSD
[ -f 'MSD_AVG' ]              && rm  MSD_AVG 
[ -f 'RDF' ]                  && rm  RDF 
[ -f 'TCF' ]                  && rm  TCF
[ -f 'TCF_AVG' ]              && rm  TCF_AVG
[ -f 'UNCHANGED_CHEMISTRY' ]  && rm  UNCHANGED_CHEMISTRY
[ -f 'RES_TIMES' ]            && rm  RES_TIMES
