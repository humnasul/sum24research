// edited by Humna Sultan for use with OC-SMART data

#ifndef _L2_FLAGS_H
#define _L2_FLAGS_H

#define NFLAGS              7
 
#define VALID               0
#define L1_REFLECTANCE      1
#define VIEW_OUT_OF_RANGE   4
#define LAND               16
#define CLOUD              64
#define LRC_OOS           256
#define LRC_NEG          1024

static char *l2_flag_lname[NFLAGS] = {"VALID",
                                      "L1_REFLECTANCE",
                                      "VIEW_OUT_OF_RANGE",
                                      "LAND",
                                      "CLOUD",
                                      "LRC_OOS",
                                      "LRC_NEG"};

static char *l2_flag_sname[NFLAGS] = {"f01_name", 
                                      "f02_name",
                                      "f03_name",
                                      "f04_name",
                                      "f05_name",
                                      "f06_name",
                                      "f07_name"};

#endif
