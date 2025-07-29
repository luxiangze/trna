#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
示例脚本：演示如何使用通用日志模块

这个脚本展示了如何在其他脚本中使用 logger_utils.py 模块
"""

import time
from logger_utils import setup_logger, get_logger, set_log_level, LOG_LEVELS


def example_function_1():
    """示例函数1：展示基本日志使用"""
    logger = get_logger('example_script')
    
    logger.info("开始执行示例函数1")
    logger.debug("这是一条调试信息（默认不显示）")
    
    # 模拟一些处理过程
    for i in range(3):
        logger.info(f"处理步骤 {i+1}/3")
        time.sleep(0.5)
    
    logger.info("示例函数1执行完成")


def example_function_2():
    """示例函数2：展示警告和错误日志"""
    logger = get_logger('example_script')
    
    logger.info("开始执行示例函数2")
    
    # 模拟警告情况
    logger.warning("这是一条警告信息：某些数据可能不完整")
    
    # 模拟错误处理
    try:
        # 故意制造一个错误
        result = 10 / 0
    except ZeroDivisionError as e:
        logger.error(f"捕获到错误: {e}")
        logger.info("错误已处理，继续执行")
    
    logger.info("示例函数2执行完成")


def main():
    """主函数"""
    # 初始化日志系统
    logger = setup_logger(__file__)
    
    logger.info("=" * 50)
    logger.info("开始执行示例脚本")
    logger.info("=" * 50)
    
    # 执行示例函数
    example_function_1()
    
    logger.info("-" * 30)
    
    example_function_2()
    
    logger.info("-" * 30)
    
    # 演示动态调整日志级别
    logger.info("演示动态调整日志级别到DEBUG")
    set_log_level(logger, LOG_LEVELS['DEBUG'])
    
    logger.debug("现在可以看到调试信息了")
    logger.info("这是调整级别后的信息")
    
    # 恢复日志级别
    logger.info("恢复日志级别到INFO")
    set_log_level(logger, LOG_LEVELS['INFO'])
    
    logger.debug("这条调试信息又看不到了")
    logger.info("脚本执行完成！")


if __name__ == '__main__':
    main()
