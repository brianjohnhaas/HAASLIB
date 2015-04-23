#!/usr/local/bin/python

import sys
import string
import Sybase


def connect_to_database (server, user, password, db):
    dbproc = Sybase.connect (server, user, password, db, auto_commit = 1 )
    return (dbproc)


def do_sql_2D (dbproc, query):
    c = dbproc.cursor()
    c.execute (query)
    return (c.fetchall())


def execute_query (dbproc, query):
    c = dbproc.cursor()
    c.execute(query)
    return None



