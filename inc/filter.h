/*
 * Copyright (c) 2020, RudyLo <luhuadong@163.com>
 *
 * SPDX-License-Identifier: MIT License
 *
 * Change Logs:
 * Date           Author       Notes
 * 2020-05-28     luhuadong    the first version
 */

struct limit_fhandle
{
    int  value;
    uint limit;
};
typedef struct limit_fhandle *limit_fhandle_t;