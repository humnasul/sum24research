#ifndef _L2_FLAGS_H
#define _L2_FLAGS_H

#define NFLAGS             32
 
#define ATMFAIL             1
#define LAND                2
#define BADANC              4
#define HIGLINT             8
#define HILT               16
#define HISENZ             32
#define COASTZ             64
#define NEGLW             128
#define STRAYLIGHT        256
#define CLDICE            512
#define COCCOLITH        1024
#define TURBIDW          2048
#define HISOLZEN         4096
#define HITAU            8192
#define LOWLW           16384
#define CHLFAIL         32768
#define NAVWARN         65536
#define ABSAER         131072
#define CLDSHDSTL      262144
#define MAXAERITER     524288
#define MODGLINT      1048576
#define CHLWARN       2097152
#define ATMWARN       4194304
#define DARKPIXEL     8388608
#define SEAICE       16777216
#define NAVFAIL      33554432
#define FILTER       67108864
#define SEAICE_ANA  134217728
#define SSTFAIL     268435456
#define NIR_SWITCH  536870912
#define FOG        1073741824
#define OCEAN      2147483648 

static char *l2_flag_lname[NFLAGS] = {"ATMFAIL", 
                                      "LAND",
                                      "BADANC",
                                      "HIGLINT",
                                      "HILT",
                                      "HISATZEN",
                                      "COASTZ",
                                      "NEGLW",
                                      "STRAYLIGHT",
                                      "CLDICE",
                                      "COCCOLITH",
                                      "TURBIDW",
                                      "HISOLZEN",
                                      "HITAU",
                                      "LOWLW",
                                      "CHLFAIL",
                                      "NAVWARN",
                                      "ABSAER",
                                      "CLDSHDSTL",
                                      "MAXAERITER",
                                      "MODGLINT",
                                      "CHLWARN",
                                      "ATMWARN",
                                      "DARKPIXEL",
                                      "SEAICE",
                                      "NAVFAIL",
                                      "FILTER",
                                      "SEAICE_ANA",
                                      "SSTFAIL",
                                      "NIR_SWITCH",
                                      "FOG",
                                      "OCEAN"};

static char *l2_flag_sname[NFLAGS] = {"f01_name", 
                                      "f02_name",
                                      "f03_name",
                                      "f04_name",
                                      "f05_name",
                                      "f06_name",
                                      "f07_name",
                                      "f08_name",
                                      "f09_name",
                                      "f10_name",
                                      "f11_name",
                                      "f12_name",
                                      "f13_name",
                                      "f14_name",
                                      "f15_name",
                                      "f16_name",
                                      "f17_name",
                                      "f18_name",
                                      "f19_name",
                                      "f20_name",
                                      "f21_name",
                                      "f22_name",
                                      "f23_name",
                                      "f24_name",
                                      "f25_name",
                                      "f26_name",
                                      "f27_name",
                                      "f28_name",
                                      "f29_name",
                                      "f30_name",
                                      "f31_name",
                                      "f32_name"};

#endif
