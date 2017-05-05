# This script takes as input the identifiers of a par set
# as well as a couple of target outcome values to change
# and then prints the insert statement for the new parameter sets
# Example call: python par_set_cloner.py --biol-prot=biol --id=0 --version=v0.4 --ps=1 --target-ps=new_parameter_set_name

import optparse, MySQLdb

db = MySQLdb.connect(user = "root", passwd = 'TCADgBq7ShmpYfN3', db = 'assay_calib_sims', 
        host = 'localhost')
con = db.cursor()

usage = """%prog [options]
It is often required to make a copy of a parameter set with only one or two very small changes.
 
This script takes as input the identifiers for a parameter set and prints to stdout the SQL commands to make a copy of that parameter set.
 
Redirect the output to a text file, modify manually and copy and paste into  MySQL shell to make slightly different copies of existing sets"""

parser = optparse.OptionParser(usage=usage)
parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False, help="Produce verbose output on stderr; Default = off")
parser.add_option("--biol-prot", action="store", dest="biol_prot", default=None, help="Cloning a Biology (biol) or Protocol (prot)?")
parser.add_option("--id", action="store", dest="id", default=None, help="Id of biology or protocol to clone")
parser.add_option("--version", action="store", dest="version", default=None, help="Version of biology or protocol to clone")
parser.add_option("--ps", action="store", dest="ps", default=None, help="Name of the parameter set to clone")
parser.add_option("--target-ps", action="store", dest="target_ps", default=None, help="Name of the new parameter set to be created")

(options, args) = parser.parse_args()

table_name = "biology_parameter_sets" if options.biol_prot == "biol" else "protocol_parameter_sets"
table_id = "biol_id" if options.biol_prot == "biol" else "prot_id"
table_ver = "biol_ver" if options.biol_prot == "biol" else "prot_ver"

query = """SELECT *
FROM {table_name}
WHERE {table_id} = {id}
AND {table_ver} = "{ver}"
AND ps = "{ps}"
""".format(table_name = table_name,
        table_id = table_id,
        table_ver = table_ver,
        id = options.id,
        ver = options.version,
        ps = options.ps)

if options.verbose:
    print query
con.execute(query)
x = con.fetchall()

print

for i in x:
    xid = i[0]
    xver = i[1]
    xps = i[2] if options.target_ps == None else options.target_ps
    xvariable = i[3]
    xtype = i[4]
    xvalue = i[5]
    xposition = i[6]
    query = """INSERT INTO {table_name} VALUES ({id}, "{ver}", "{ps}", "{variable}", "{type}", "{value}", {position});""".format(table_name = table_name,
            id = xid, 
            ver = xver, 
            ps = xps, 
            variable = xvariable, 
            type = xtype, 
            value = xvalue,
            position = xposition)
    print query
