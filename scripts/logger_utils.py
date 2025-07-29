#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
通用日志工具模块

提供统一的日志配置功能，支持多脚本使用。
日志文件自动保存到 logs/<脚本名>/日期时间.log

使用方法：
    from logger_utils import setup_logger
    
    # 在脚本开始处调用
    logger = setup_logger(__file__)
    
    # 使用日志
    logger.info("这是一条信息")
    logger.warning("这是一条警告")
    logger.error("这是一条错误")
"""

import os
import logging
from datetime import datetime
from pathlib import Path


def setup_logger(script_path, log_level=logging.INFO, console_output=True):
    """
    设置日志配置，为指定脚本创建专用的日志系统
    
    Args:
        script_path (str): 脚本文件路径，通常传入 __file__
        log_level (int): 日志级别，默认为 logging.INFO
        console_output (bool): 是否同时输出到控制台，默认为 True
    
    Returns:
        logging.Logger: 配置好的日志器对象
    """
    # 获取脚本名称（不含扩展名）
    script_name = Path(script_path).stem
    
    # 创建日志目录
    log_dir = os.path.join(os.getcwd(), 'logs', script_name)
    os.makedirs(log_dir, exist_ok=True)
    
    # 生成日志文件名
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_file = os.path.join(log_dir, f'{timestamp}.log')
    
    # 创建日志器
    logger = logging.getLogger(script_name)
    
    # 清除已有的处理器（避免重复配置）
    if logger.handlers:
        logger.handlers.clear()
    
    # 设置日志级别
    logger.setLevel(log_level)
    
    # 创建格式器
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # 创建文件处理器
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    # 创建控制台处理器（可选）
    if console_output:
        console_handler = logging.StreamHandler()
        console_handler.setLevel(log_level)
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)
    
    # 防止日志传播到根日志器
    logger.propagate = False
    
    # 记录日志系统启动信息
    logger.info(f"日志系统已启动 - 脚本: {script_name}")
    logger.info(f"日志文件: {log_file}")
    logger.info(f"日志级别: {logging.getLevelName(log_level)}")
    logger.info(f"控制台输出: {'启用' if console_output else '禁用'}")
    
    return logger


def get_logger(script_name):
    """
    获取已配置的日志器
    
    Args:
        script_name (str): 脚本名称
    
    Returns:
        logging.Logger: 日志器对象，如果未配置则返回基础日志器
    """
    return logging.getLogger(script_name)


def set_log_level(logger, level):
    """
    动态设置日志级别
    
    Args:
        logger (logging.Logger): 日志器对象
        level (int): 新的日志级别
    """
    logger.setLevel(level)
    for handler in logger.handlers:
        handler.setLevel(level)
    logger.info(f"日志级别已更改为: {logging.getLevelName(level)}")


def add_file_handler(logger, log_file, level=logging.INFO):
    """
    为现有日志器添加额外的文件处理器
    
    Args:
        logger (logging.Logger): 日志器对象
        log_file (str): 日志文件路径
        level (int): 日志级别
    """
    formatter = logging.Formatter(
        '%(asctime)s - %(levelname)s - %(funcName)s:%(lineno)d - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    file_handler = logging.FileHandler(log_file, encoding='utf-8')
    file_handler.setLevel(level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    logger.info(f"已添加额外日志文件: {log_file}")


# 便捷的日志级别常量
LOG_LEVELS = {
    'DEBUG': logging.DEBUG,
    'INFO': logging.INFO,
    'WARNING': logging.WARNING,
    'ERROR': logging.ERROR,
    'CRITICAL': logging.CRITICAL
}


if __name__ == '__main__':
    # 测试代码
    logger = setup_logger(__file__)
    
    logger.debug("这是一条调试信息")
    logger.info("这是一条普通信息")
    logger.warning("这是一条警告信息")
    logger.error("这是一条错误信息")
    logger.critical("这是一条严重错误信息")
    
    print(f"日志文件保存在: logs/{Path(__file__).stem}/")
