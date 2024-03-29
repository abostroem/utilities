import os
from astropy.io import ascii as asc
import pymysql

def database_config(config_file):
    config = asc.read(config_file, data_start=0)
    config_dict = {}
    for keyword, arg in config:
        config_dict[keyword]=arg
    return config_dict

def get_cursor(config_filename='configure', config_dir='/Users/bostroem/lco'):
    """
    Connect to a database specified in <config_dir>/<config_filename>
    e.g. local supernova: use defaults; dark supernova: config_filename='configure_dark'
    """
    config_dict = database_config(os.path.join(config_dir, config_filename))
    if 'port' in config_dict.keys():
        db = pymysql.connect(user=config_dict['mysqluser'], 
                        host=config_dict['hostname'], 
                        password=config_dict['mysqlpasswd'],
                        database=config_dict['database'],
                    port=int(config_dict['port']),
                    cursorclass=pymysql.cursors.DictCursor)
    else:
        db = pymysql.connect(user=config_dict['mysqluser'], 
                        host=config_dict['hostname'], 
                        password=config_dict['mysqlpasswd'],
                        database=config_dict['database'],
                        cursorclass=pymysql.cursors.DictCursor)
    cursor = db.cursor()
    return db, cursor