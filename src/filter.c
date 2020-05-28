/*
 * Copyright (c) 2020, RudyLo <luhuadong@163.com>
 *
 * SPDX-License-Identifier: MIT License
 *
 * Change Logs:
 * Date           Author       Notes
 * 2020-05-28     luhuadong    the first version
 */

#include <rtthread.h>
#include "filter.h"

/*
 * 名称：限幅滤波法（又称程序判断滤波法）
 * 方法：
 *  根据经验判断，确定两次采样允许的最大偏差值（设为A），
 *  每次检测到新值时判断：
 *  如果本次值与上次值之差<=A，则本次值有效，
 *  如果本次值与上次值之差>A，则本次值无效，放弃本次值，用上次值代替本次值。
 * 优点：
 *  能有效克服因偶然因素引起的脉冲干扰。
 * 缺点：
 *  无法抑制那种周期性的干扰。
 *  平滑度差。
 * 说明：
 *  限幅滤波法主要用于处理变化较为缓慢的数据，如温度、物体的位置等。
 *  使用时，关键要选取合适的门限制A。通常这可由经验数据获得，必要时可通过实验得到。
*/
int limit_filter(limit_fhandle_t handle, int value)
{
    RT_ASSERT(handle);

    int delta = abs(value - handle->value);

    if (delta < handle->limit)
        handle->value = value;
    
    return handle->value;
}

int median_filter()
{

}

int mean_filter()
{

}

int weighted_mean_filter()
{

}

int move_mean_filter()
{

}

int low_pass_filter()
{

}

int high_pass_filter()
{
    
}