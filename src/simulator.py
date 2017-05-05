#!/usr/bin/python

# This code will simulate patient's biomarker values over time.

# The values of the biomarkers are controlled by a biology (also referred to as
# a functional form) and the timing of
# the measurements are controlled a protocol. Hence, in order to specify a simulation 
# you must specify a biology and a protocol.

# Details of how to use the code can be found by running python simulator.py -h
# or by reading the parser.add_option commands in this script.
# example call: python simulator.py -c sims -p "prot3.1m_1m_19_10" -b "biology2.noise_015_bed" -n 10

# The start of a numerical integration based approach to computing exact
# solutions for the mdr of certain biologies are in this code. This was
# scrapped and moved to R due to more familiarity with R's integration code.
# However, there is a trick which can be used to use this simulation code to
# estimate the exact solutions. By running a protocol with perfect daily
# followup and a very large number of patients, you can simulate enough data
# that you can almost exactly compute the mdr of the biology - however the
# numerical integration approach is much better.

# There used to be a WebUI that interfaced with this script. It used an object
# orientated approach and hence there are several classes and methods in this
# script that still contains code related to these functions. I did not have
# time to remove them.

# The code is organized using folding with the fold markers {{{ and }}}.
# If you use an editor like vim you can set foldmethod=marker and then the
# structure of the code will be clearer. The basic outline is:
# imports: load python modules
# Database Connection 
# Commanline Argument Parsing
# Class Declarations - much of this is focussed on the webui.
# Initialise database tables
# Specify the Function Forms (Biologies)
# Specify Protocols 
# Functions to save data in the database 
# The big function that simulates the data: simulateCohorts
# A function related to the computation of exact solutions: exactSolutions
# Prepare the biologies and protocols from the commandline arguments

# There are more comments in the sections on the biologies and protocols that
# contain instructions on how to construct new ones.

# imports: load python modules {{{
import datetime
import random
import math
#import utilities as util
#import simulator_cfg as cfg
from UserDict import UserDict
import MySQLdb
import scipy.stats
from optparse import OptionParser
import numpy
#from cspy.utility.html import TextSetting
import numpy as np
import time
#import logging
#logging.basicConfig(level=logging.debug, filename = 'simulator.log') #}}}

# Database Connection {{{
db = MySQLdb.connect(user = "root", passwd = 'TCADgBq7ShmpYfN3', db = 'assay_calib_sims', host = 'localhost')
con = db.cursor() #}}}

# Commanline Argument Parsing {{{
parser = OptionParser()
parser.add_option("-c", "--command", dest="command", help="What action to perform: 'sims' for simulating a cohort's biomarker value or 'exact' for computing an exact solution.")
parser.add_option("-p", "--protocols", dest="reqProtocols", help="What protocols and parameter sets to use. Notation: The name of the protocol followed by a dot '.' followed by a list of the parameter sets seperated by semicolons ';'. For example: protx.1m_1w_1w;1m_1m_1w. Mutliple protocols can be specified for a single run, but this is generally not a good idea. The name of the protocol must be in the protocolD dictionary defined in the script. The names of the parameter sets must be in the MySQL database that the script connects to in the protocol_parameter_sets table.")
parser.add_option("-b", "--biologies", dest="reqBiologies", help="What biologies and parameter sets to use. See the description of protocols for details of the notation. Biologies are defined in the biologiesD dictionary and the biology_parameter_sets table.")
parser.add_option("-n", "--ncohorts", dest="ncohorts_input", help="How many cohorts to simulate for each scenario")
parser.add_option("--do-not-round", dest="roundVisitDates", default = True, action = "store_false", help="Should the visit dates be rounded?")
parser.add_option("--do-not-restrict-bmv", dest="restrictBMV", default = True, action = "store_false", help="Should the BMVs be resticted to be greater than zero?")
parser.add_option("--threshold", dest="threshold", default = 0, action = "store", type = "float", help="The threshold for the exact solution. WARNING: computation of exact methods are buggy - rather use R")
parser.add_option("--bigT", dest="bigT", default = 365, action = "store", type = "float", help="The bigT for the exact solution. WARNING: computation of exact methods are buggy - rather use R")
parser.add_option("--integration-method", dest="integration_method", default = "quad", action = "store", type = "string", help="Which integration method should be used for the computation of the exact method? quad / monte carlo. WARNING: computation of exact methods are buggy - rather use R")
(options, args) = parser.parse_args()
#}}}

# Class Declarations {{{

class ReportOption(object): #{{{@never_cache

    def __init__(self,
            name,
            title = False,
            html = False,
            default = False,
            altvalues = False,
            itype = False,
            args = False,
            help_text = False,
            required = False,
            validators=[]):
                
        self.name = name
        self.title = title
        self.default = default
        self.altvalues = altvalues
        if html:
            self.html = html
        self.itype=itype
        if args:
            self.value = self.get_value(args)
        else:
            self.value = False
        self.help_text = help_text
        self.error = False
        self.required = required
        self.validators = validators

    def get_value(self, args):
        value = args.get(self.name, self.default)
        if value == "":
            value = self.default
        return value
    
    def get_all(self):
        return [self.default] + self.altvalues

    def get_default(self, args):
        return self.default

    def validate(self, args):
        value = self.get_value(args)
        self.error=""
        for validator in self.validators:
            try:
                value, error = validator(value)
            except ValidationError as e:
                print "ValidationError caught:", e.messages[0]
                value = value
                error = e.messages[0]
            if error:
                self.error = error
                break
        self.value = value
        #print "validation value:", value
        return value

    def html_title(self):
        if not self.title:
            return ""
        return '''<small><b>%s</b></small><br/>\n''' % self.title

    def html_help_text(self):
        if self.help_text:
            return '''<div class="help_text">%s</div>\n''' % self.help_text
        else:
            return ""

    def html_errors(self):
        html = ""
        if self.error:
            html+='<div class="error_message">%s</div>\n' %self.error
        return html

    def html(self):
        return self.html_title() + self.html_help_text() + self.html_errors() + self.html_input() + "<br/>"
#}}}

class TextSetting(ReportOption):#{{{
    '''Class to display a text input element in an html page and parse the response from the user.'''

    def __init__(self, name, default="", *args, **kwargs):
        '''Compulsary args:
        name - string - python safe and html safe name that will serve as the name of the input element and the key in the args dictionary.
        title - string - User friendly string that is displayed above the input element, telling the user what setting he is changing.
        Optional args:
        default - string - Value that will be entered into text-field the first time it is rendered.'''
        ReportOption.__init__(self, name=name, default=default, *args, **kwargs)
        self.itype = "textfield"
        
    def html_input(self):
        '''returns html string for the input element.'''
        if self.value:
            value = self.value
        else:
            value = self.default
        html = '''<input type="text" name = "%s" value = "%s"/><br/>\n''' %(self.name, value)
        return html #}}}

class AnnotatedTextSetting(TextSetting): #{{{
    def __init__(self, db_var_name, db_var_type = "double", *args, **kwargs):
        TextSetting.__init__(self, *args, **kwargs)
        self.db_var_name = db_var_name 
        self.db_var_type = db_var_type#}}}

class odict(UserDict): #{{{
    """
    Use this class to pass the parameters
    The key based storage system makes unpacking the parameters explicit
    The fact that it is ordered makes writing the values to files easier
    """
    def __init__(self, dict = None):
        self._keys = []
        UserDict.__init__(self, dict)

    def __delitem__(self, key):
        UserDict.__delitem__(self, key)
        self._keys.remove(key)

    def __setitem__(self, key, item):
        UserDict.__setitem__(self, key, item)
        if key not in self._keys: self._keys.append(key)

    def clear(self):
        UserDict.clear(self)
        self._keys = []

    def copy(self):
        dict = UserDict.copy(self)
        dict._keys = self._keys[:]
        return dict

    def items(self):
        return zip(self._keys, self.values())

    def keys(self):
        return self._keys

    def popitem(self):
        try:
            key = self._keys[-1]
        except IndexError:
            raise KeyError('dictionary is empty')

        val = self[key]
        del self[key]

        return (key, val)

    def setdefault(self, key, failobj = None):
        UserDict.setdefault(self, key, failobj)
        if key not in self._keys: self._keys.append(key)

    def update(self, dict):
        UserDict.update(self, dict)
        for key in dict.keys():
            if key not in self._keys: self._keys.append(key)

    def values(self):
        return map(self.get, self._keys) #}}}

class ParSet(object):#{{{

    def __init__(self, par_type, type_id, type_ver, ps, parameters):
        self.par_type = par_type
        self.type_id = type_id
        self.type_ver = type_ver
        self.ps = ps
        self.parameters = parameters

    def filters(self, prefix):
        """
        This method should generate the filters for the parameter set
        It should add a prefix to the names of the filters
        Generating the Titles are still n problem, it should
            be read from the databse? Not sure
        Note that in simulateNewDataset.py, the syntax of the name is used to
            deduce the names of the columns in the database
        """
        pass

    #}}}

class Biology(odict):#{{{

    def __init__(self, id, version, table, parameters, dict=None):#{{{
        """
        table is an object of class BiologyParameterSetTable
        """
        odict.__init__(self)
        self.id = id
        self.version = version
        self.table = table
        self.parameters = parameters
        self.refresh_par_sets()
    #}}}

    def refresh_par_sets(self): #{{{
        #print "refreshing parsets"
        self.par_set_ids = self.table.get_par_set_ids(self.id, self.version)
        self.par_sets = {}
        for par_set_id in self.par_set_ids:
            self.par_sets[par_set_id] = self.table.load_set_as_object(biol_id = self.id, 
                    biol_ver = self.version, 
                    ps = par_set_id)
        #print "new parsets:"
        #print self.par_sets
        #}}}

    def par_set_formatted(self, par_set_id): #{{{
        html = ""
        for i in range(len(self.parameters)):
            html += '''<small><b>%s</b></small><br/>
            %s<br/>\n''' %(self.parameters[i].title, self.par_sets[par_set_id].parameters[i]['value'])
            html += '''<div  style="margin: 0px; padding:0px; clear: both; height: 16px">&nbsp;</div>\n'''
        return html
    #}}}

#}}}

class Protocol(odict):#{{{

    def __init__(self, id, version, table, parameters, dict=None):#{{{
        """
        table is an object of type ProtocolParameterSetTable
        """
        odict.__init__(self)
        self.id = id
        self.version = version
        self.table = table
        self.parameters = parameters
        self.refresh_par_sets()
        #}}}

    def refresh_par_sets(self): #{{{
        #print "refreshing parsets"
        self.par_set_ids = self.table.get_par_set_ids(self.id, self.version)
        self.par_sets = {}
        for par_set_id in self.par_set_ids:
            self.par_sets[par_set_id] = self.table.load_set_as_object(prot_id = self.id, 
                    prot_ver = self.version, 
                    ps = par_set_id)
        #print "new parsets:"
        #print self.par_sets
        #}}}

    def par_set_formatted(self, par_set_id): #{{{
        html = ""
        for i in range(len(self.parameters)):
            html += '''<small><b>%s</b></small><br/>
            %s<br/>\n''' %(self.parameters[i].title, self.par_sets[par_set_id].parameters[i]['value'])
            html += '''<div  style="margin: 0px; padding:0px; clear: both; height: 16px">&nbsp;</div>\n'''
        return html
    #}}}
    #}}}
#}}}  end of class definitions

# Initialise database tables {{{

class BiologyParameterSetTable(object):#{{{

    def __init__(self):
        self.tableName = 'biology_parameter_sets'
        pass

    def create_table(self):#{{{
        con.execute("DROP TABLE IF EXISTS `%s`" %(self.tableName))
        con.execute('''CREATE TABLE `%s` (
          `biol_id` int(11) NOT NULL,
          `biol_ver` varchar(45) DEFAULT NULL,
          `ps` int(11) DEFAULT NULL,
          `variable` varchar(45) DEFAULT NULL,
          `type` varchar(45) DEFAULT NULL,
          `value` varchar(45) DEFAULT NULL,
          `position` INTEGER DEFAULT NULL
          ) ENGINE=InnoDB''' %(self.tableName) )
        db.commit()
        return True
        #}}}

    def save_set(self, new_set):#{{{
        position = 1
        for par in new_set["parameters"]:
            query = '''INSERT INTO %s (biol_id, biol_ver, ps, variable, type, value, position)
            VALUES ( '%s', '%s', '%s', '%s', '%s', '%s', '%s' )''' %(self.tableName,
                new_set["biol_id"],
                new_set["biol_ver"],
                new_set["ps"],
                par["variable"],
                par["type"],
                par["value"],
                str(position))
            #print query
            con.execute(query)
            position += 1
        db.commit()
        return True
        #}}}

    def get_par_set_ids(self, biol_id = 'False == 0 is True - darn', biol_ver = 'False == 0 is True - darn'):#{{{
        # The default values assigned might hint at a bug that was tricky to catch
        query = '''SELECT distinct(ps) FROM %s a
            INNER JOIN (select biol_id pid, max(biol_ver) pvr 
            from %s 
            group by biol_id) b
            on a.biol_id = b.pid and a.biol_ver = b.pvr\n'''%(self.tableName, self.tableName)
        if biol_id != 'False == 0 is True - darn' and biol_ver != 'False == 0 is True - darn':
            query += '''WHERE biol_id = "%s"
            AND biol_ver = "%s"\n''' %(biol_id, biol_ver)
        query += '''ORDER BY ps''' 
        con.execute(query)
        data = con.fetchall()
        db.commit()
        return [ r[0] for r in data ]
    #}}}

    def load_set(self, biol_id, biol_ver, ps):#{{{
        query = '''SELECT variable, type, value, position from %s
        WHERE biol_id = "%s"
        AND biol_ver = "%s"
        AND ps = '%s'
        ORDER BY position''' %(self.tableName, biol_id, biol_ver, ps)
        con.execute(query)
        data = con.fetchall()
        par_set = {"biol_id": biol_id,
            "biol_ver": biol_ver,
            "ps": ps,
            "parameters": [    # Begin list comprehension
                    {"variable": par[0],
                    "type": par[1],
                    "value": par[2]}
                for par in data
                ]
            }
        x = par_set
        combo = odict()
        combo['bio_param_set'] = x['ps']
        for i in x['parameters']:
            if i['type'] == "double":
                combo[i['variable']] = float(i['value'])
            else:
                combo[i['variable']] = str(i['value'])
        return [par_set, combo]
        #}}} 
    
    def find_similar_set(self, new_set):#{{{
        """
        Given a new set of parameters,
            This function will find all ps_id's
            extract each of those ps's
            make sure that the new_set is unique (w.r.t. parameter values and not par_set_id)
            and return -1 if unique
            OR
            return the par_set_id (ps column) of the param set that is identical to the currrent one
        """
        query = '''SELECT DISTINCT ps from %s
        WHERE biol_id = '%s'
        AND biol_ver = '%s'
        ''' %(self.tableName, new_set['biol_id'], new_set['biol_ver'])
        con.execute(query)
        data = con.fetchall()
        for i in data:
            the_same = True
            x = self.load_set(new_set['biol_id'], new_set['biol_ver'], i[0])[0]
            for j in x['parameters']:
                for l in new_set['parameters']:
                    if j['variable'] == l['variable']:
#                        print "match"
#                        print j['variable']
#                        print j['value'], l['value']
                        if j['value'] != l['value']:
                            the_same = False
#                        print '-'*80
            if the_same:
#                print "similar set found"
#                print i
                return i[0]
            else:
                return -1

        #}}}

    def load_set_as_object(self, biol_id, biol_ver, ps):#{{{
        as_dict = self.load_set(biol_id, biol_ver, ps)[0]
        parameters = as_dict["parameters"]
        new_obj = ParSet(par_type="biol", type_id = biol_id, type_ver = biol_ver, ps = ps, parameters = parameters)
        return new_obj
        #}}}

bioObj = BiologyParameterSetTable()
    #}}}

class ProtocolParameterSetTable(object):#{{{

    def __init__(self):
        self.tableName = 'protocol_parameter_sets'
        pass

    def create_table(self):#{{{
        con.execute("DROP TABLE IF EXISTS `%s`" %(self.tableName))
        con.execute('''CREATE TABLE `%s` (
          `prot_id` int(11) NOT NULL,
          `prot_ver` varchar(45) DEFAULT NULL,
          `ps` int(11) DEFAULT NULL,
          `variable` varchar(45) DEFAULT NULL,
          `type` varchar(45) DEFAULT NULL,
          `value` varchar(45) DEFAULT NULL,
          `position` INTEGER DEFAULT NULL
          ) ENGINE=InnoDB''' %(self.tableName) )
        db.commit()
        return True
        #}}}

    def save_set(self, new_set):#{{{
        """
        new_set is a dictionary containing the values of the identifiers
        AND a dictionary of the parameters and their values
        """
        position = 1
        for par in new_set["parameters"]:
            query = '''INSERT INTO %s (prot_id, prot_ver, ps, variable, type, value, position)
            VALUES ( '%s', '%s', '%s', '%s', '%s', '%s', '%s' )''' %(self.tableName,
                new_set["prot_id"],
                new_set["prot_ver"],
                new_set["ps"],
                par["variable"],
                par["type"],
                par["value"],
                str(position))
            #print query
            con.execute(query)
            position += 1
        db.commit()
        return True
        #}}}

    def get_par_set_ids(self, prot_id='False == 0 is True - darn', prot_ver='False == 0 is True - darn'):#{{{
        query = '''SELECT distinct(ps) FROM %s a 
            INNER JOIN (select prot_id pid, max(prot_ver) pvr 
            from %s 
            group by prot_id) b
            on a.prot_id = b.pid and a.prot_ver = b.pvr\n'''%(self.tableName, self.tableName)
        if prot_id != 'False == 0 is True - darn' and prot_ver != 'False == 0 is True - darn':
            query += '''WHERE prot_id = "%s"
            AND prot_ver = "%s"\n''' %(prot_id, prot_ver)
        query += '''ORDER BY ps''' 
        con.execute(query)
        data = con.fetchall()
        db.commit()
        return [ r[0] for r in data ]
    #}}}

    def load_set(self, prot_id, prot_ver, ps):#{{{
        query = '''SELECT variable, type, value, position from %s
        WHERE prot_id = "%s"
        AND prot_ver = "%s"
        AND ps = "%s"
        ORDER BY position''' %(self.tableName, prot_id, prot_ver, ps)
        con.execute(query)
        data = con.fetchall()
        par_set = {"prot_id": prot_id,
            "prot_ver": prot_ver,
            "ps": ps,
            "parameters": [
                {"variable": par[0],
                "type": par[1],
                "value": par[2]}
                for par in data
                ]
            }
        x = par_set
        combo = odict()
        combo['prot_param_set'] = x['ps']
        for i in x['parameters']:
            if i['type'] == "Integer":
                combo[i['variable']] = int(i['value'])
            else:
                combo[i['variable']] = float(i['value'])
        return [par_set, combo]
                
        #}}}

    def find_similar_set(self, new_set):#{{{
        """
        Given a new set of parameters,
            This function will find all ps_id's
            extract each of those ps's
            make sure that the new_set is unique (w.r.t. parameter values and not par_set_id)
            and return -1 if unique
            OR
            return the par_set_id (ps column) of the param set that is identical to the currrent one
        """
        query = '''SELECT DISTINCT ps from %s
        WHERE prot_id = '%s'
        AND prot_ver = '%s'
        ''' %(self.tableName, new_set['prot_id'], new_set['prot_ver'])
        con.execute(query)
        data = con.fetchall()
        for i in data:
            the_same = True
            x = self.load_set(new_set['prot_id'], new_set['prot_ver'], i[0])[0]
            for j in x['parameters']:
                for l in new_set['parameters']:
                    if j['variable'] == l['variable']:
#                        print "match"
#                        print j['variable']
#                        print j['value'], l['value']
                        if j['value'] != l['value']:
                            the_same = False
#                        print '-'*80
            if the_same:
#                print "similar set found"
#                print i
                return i[0]
            else:
                return -1

    def load_set_as_object(self, prot_id, prot_ver, ps):#{{{
        as_dict = self.load_set(prot_id, prot_ver, ps)[0]
        parameters = as_dict["parameters"]
        new_obj = ParSet(par_type="prot", type_id = prot_id, type_ver = prot_ver, ps = ps, parameters = parameters)
        return new_obj
        #}}}

protObj = ProtocolParameterSetTable()
 #}}}

# Phils parameter extractor
# dunno if this should go into the class...?

def extractPSbio(biol):#{{{
    id = biol['biol_id']
    version = biol['version']
    pss = []
    for i in range(1000):
        try:
            pss.append(biol['cohort_pars_gen'](1, 1, i))
        except:
            pass
    for i in pss:
        new_set = odict()
        new_set['biol_id'] = id
        new_set['biol_ver'] = version
        new_set['ps'] = i['biol_param_set']
        new_set['parameters'] = []
        for j in i.keys():
            if j != "biol_param_set":
                par = odict()
                par['variable'] = j
                par['value'] = i[j]
                par['type'] = "Double"
                new_set['parameters'].append(par)
        #print new_set
        bioObj.save_set(new_set)
# sample usage
# extractPSbio(biologiesD['biology1'])
# However, this will only work if you return the combo1, combo2 functionality to how it was  #}}} 

def extractPSpro(prot):#{{{
    id = prot['prot_pars_gen']()['prot_id']
    version = prot['prot_pars_gen']()['prot_version']
    pss = []
    for i in range(1000):
        try:
            pss.append(prot['cohort_pars_gen'](1, i))
        except:
            pass
    for i in pss:
        new_set = odict()
        new_set['prot_id'] = id
        new_set['prot_ver'] = version
        new_set['ps'] = i['prot_param_set']
        new_set['parameters'] = []
        for j in i.keys():
            if j != "prot_param_set":
                par = odict()
                par['variable'] = j
                par['value'] = i[j]
                par['type'] = "Double"
                new_set['parameters'].append(par)
        #print new_set
        protObj.save_set(new_set)   

    #}}}

#}}}

# }}}

# Specify the Function Forms (Biologies){{{

# Description {{{
# The function forms specify the biology of how the biomarker behaves
# The biology varies from subject to subject, but also from population to population
# A function form is described by contructing a function form family by:
#   a) Specifying the function that describes shape of curve according to 
#       which the biomarker evolve. 
#       (ffx_bmf_fun) (see below)
#   b) Specifying the distributions from which the population specific
#       parameters for the family are sampled
#       (ffx_cohort_pars_gen) (see below)
#   c) Specifying the distributions from which the subject specific 
#       parameters for the family are sampled
#       (ffx_sub_pars_gen) (see below)
# This is accompished by writing the following three functions:
#   ffx_bmf_fun
#       Given all information, it simulates the value of the biomarker
#       i.e. it generates the visit specific parameters (which is only the value of the biomarker)
#   ffx_cohort_pars_gen
#       This function simulates the values of the parameters that will be the same for
#       all the subjects in the cohort. 
#       i.e. it generates the population specific parameters
#   ffx_sub_pars_gen
#       This function simulates the values of the parameters that will be the same for
#       each subject. For example, the subject's seroconversion date
#       i.e. it generates the subject specific parameters

# A function form family must be combined with a protocol inorder for a simulation
#   to be completely specified }}}

# Specify ff1 {{{ root
def ff1_bmf_fun(ff_cohort_pars, 
        ff_sub_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        visits): #{{{

    t = visits[-1][1]
    alpha = ff_sub_pars['alpha']
    theta = ff_sub_pars['theta']
    seroconversion_date = prot_sub_pars['seroconversion_date'] 
    sigma = ff_cohort_pars['sigma']
    
    error_term = random.gauss(0, sigma)
    Z = 8 * pow(max(0, (t - seroconversion_date - alpha)),theta)
    
    if options.restrictBMV:
        return max(0, Z + error_term) 
    else:
        return Z + error_term
    #}}}

def ff1_cohort_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        param_set, 
        version,
        biol_id,
        fun = True): #   {{{

    return bioObj.load_set(biol_id, version, param_set)[1]

def ff1_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id, 
        fun = True): #{{{

    theta_l = ff_cohort_pars['theta_l']
    theta_u = ff_cohort_pars['theta_u']
    alpha_l = ff_cohort_pars['alpha_l']
    alpha_u = ff_cohort_pars['alpha_u']

    def alpha_gen():
        return random.uniform(alpha_l, alpha_u)

    def theta_gen():
        return random.uniform(theta_l, theta_u)

    results = odict()
    results['alpha'] = alpha_gen()
    results['theta'] = theta_gen()

    return results #}}}

biology1 = Biology(id=0, version="v0.6", table=bioObj, parameters = [#{{{
        AnnotatedTextSetting(name="biol0_par_theta_l", title="Theta: Lower Bound", default = 0.275, db_var_name = "theta_l"),
        AnnotatedTextSetting(name="biol0_par_theta_u", title="Theta: Upper Bound", default = 0.325, db_var_name = "theta_u"),
        AnnotatedTextSetting(name="biol0_par_alpha_l", title="Alpha: Lower Bound", default = 15, db_var_name = "alpha_l"),
        AnnotatedTextSetting(name="biol0_par_alpha_u", title="Alpha: Upper Bound", default = 25, db_var_name = "alpha_u"),
        AnnotatedTextSetting(name="biol0_par_sigma", title="Sigma", default = 4, db_var_name = "sigma")
        ]) 
biology1['bmf_fun'] = ff1_bmf_fun
biology1['cohort_pars_gen'] = ff1_cohort_pars_gen
biology1['sub_pars_gen'] = ff1_sub_pars_gen
biology1['version'] = 'v0.6'
biology1['biol_id'] = 0 #}}}

## }}}

# Specify ff2 {{{ log

def ff2_bmf_fun(ff_cohort_pars, 
        ff_sub_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        visits): #{{{

    t = visits[-1][1]
    seroconversion_date = prot_sub_pars['seroconversion_date']
    alpha = ff_sub_pars['alpha']
    theta = ff_sub_pars['theta']
    sigma = ff_cohort_pars['sigma']
    
    error_term = random.gauss(0, sigma)
    Z = max(0, math.log(max(t - seroconversion_date - alpha, 1), theta))
    
    if options.restrictBMV:
        return max(0, Z + error_term) 
    else:
        return Z + error_term
    #}}}

def ff2_cohort_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        param_set, 
        version,
        biol_id,
        fun = True): #{{{

    return bioObj.load_set(biol_id, version, param_set)[1]
#    param_sets = odict()
#
#    combo1 = odict()
#    combo1['biol_param_set'] = 1
#    combo1['theta_l'] = 2
#    combo1['theta_u'] = 3
#    combo1['alpha_l'] = 25
#    combo1['alpha_u'] = 35
#    combo1['sigma'] = 0.05
#    combo1['assay_threshold'] = 5
#
#    param_sets[1] = combo1
#
#    return param_sets[param_set] #}}}

def ff2_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id, 
        fun = True): #{{{

    theta_l = ff_cohort_pars['theta_l']
    theta_u = ff_cohort_pars['theta_u']
    alpha_l = ff_cohort_pars['alpha_l']
    alpha_u = ff_cohort_pars['alpha_u']

    def alpha_gen():
        return random.uniform(alpha_l, alpha_u)

    def theta_gen():
        return random.uniform(theta_l, theta_u)

    results = odict()
    results['alpha'] = alpha_gen()
    results['theta'] = theta_gen()

    return results #}}}

def ff2_exact(ff_cohort_pars, bigT = 600): #{{{
    """
    This function computes the exact solution for
    the function form by means of integration
    Important inputs are:
        parameters for distributions of each parameter
        The distributions themselves are hard-coded
    """
    theta_l = ff_cohort_pars['theta_l']
    theta_u = ff_cohort_pars['theta_u']
    alpha_l = ff_cohort_pars['alpha_l']
    alpha_u = ff_cohort_pars['alpha_u']
    sigma = ff_cohort_pars['sigma']
    assay_threshold = ff_cohort_pars['assay_threshold']
 
    # in the end it all boils down to the noise. Whats the probability of a normal variable 
    # being > (threshold - true biomarker value)
    # P(z > x) = 1 - scipy.stats.norm.cdf(x)

    def pR(t, threshold, alpha, theta, sigma):
        x = scipy.stats.norm.cdf(((threshold - math.log(max(1,(t - alpha)), theta))/sigma))
#        print (t, x)
        return x

#    def sigma_pdf(simga, ub=0.20, lb=0.1):
#        return 1/(float(ub) - lb)
#
#    def pR_1(sigma, threshold, alpha, theta, bigT, sigma_pdf):
#        #print (theta, alpha, sigma)
#        return sigma_pdf(sigma) * scipy.integrate.quad(pR, 0, bigT, args = (threshold, alpha, theta, sigma))[0]

    def alpha_pdf(alpha, ub=alpha_u, lb=alpha_l):
        return 1/(float(ub)-lb)

    def pR_2(alpha, threshold, theta, sigma, bigT, alpha_pdf):
        x = scipy.integrate.quad(pR, 0, bigT, args = (threshold, alpha, theta, sigma))[0]
        y = alpha_pdf(alpha)
#        print (sigma, theta, alpha, x, y)
        return x * y 

    def theta_pdf(theta, ub=theta_u, lb=theta_l):
        return 1/(float(ub) - lb)

    def pR_3(theta, threshold, sigma, bigT, alpha_pdf, theta_pdf):
        x = scipy.integrate.quad(pR_2, alpha_l, alpha_u, args = (threshold, theta, sigma, bigT, alpha_pdf))[0]
        y = theta_pdf(theta)
#        print (theta, x, y)
        return theta_pdf(theta) * x

    def pR_4(threshold, sigma, bigT, alpha_pdf, theta_pdf):
        return scipy.integrate.quad(pR_3, theta_l, theta_u, args = (threshold, sigma, bigT, alpha_pdf, theta_pdf))

    mdr = pR_4(assay_threshold,
            sigma,
            bigT,
            alpha_pdf,
            theta_pdf)

    return mdr #}}}

biology2 = Biology(id=1, version="v0.2", table = bioObj, parameters = [ #{{{
    AnnotatedTextSetting(name="biol1_par_theta_l", title="Theta: Lower Bound", default = 2, db_var_name = "theta_l"),
    AnnotatedTextSetting(name="biol1_par_theta_u", title="Theta: Upper Bound", default = 3, db_var_name = "theta_u"),
    AnnotatedTextSetting(name="biol1_par_alhpa_l", title="Alpha: Lower Bound", default = 25, db_var_name = "alpha_l"),
    AnnotatedTextSetting(name="biol1_par_alhpa_u", title="Alpha: Upper Bound", default = 35, db_var_name = "alpha_u"),
        AnnotatedTextSetting(name="biol1_par_sigma", title="Sigma", default = 0.05, db_var_name = "sigma"),
        AnnotatedTextSetting(name="biol1_par_assay_threshold", title="Assay Threshold", default = 5, db_var_name = "assay_threshold")
        ] )
biology2['bmf_fun'] = ff2_bmf_fun
biology2['cohort_pars_gen'] = ff2_cohort_pars_gen
biology2['sub_pars_gen'] = ff2_sub_pars_gen
biology2['exact'] = ff2_exact
#biology2['version'] ='v0.2'
biology2['biol_id'] = 1 #}}}

# }}}

# Specify ff3 {{{ weibull
def ff3_bmf_fun(ff_cohort_pars, 
        ff_sub_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        visits): #{{{

    t = visits[-1][1]
    alpha = ff_sub_pars['alpha']
    beta = ff_sub_pars['beta']
    height = ff_sub_pars['height']
    seroconversion_date = prot_sub_pars['seroconversion_date']
    ea = ff_cohort_pars['ea']
    eb = ff_cohort_pars['eb']
    ec = ff_cohort_pars['ec']
    ed = ff_cohort_pars['ed']

    Z = height*(1 - math.exp(-(math.exp(-beta)*(t - seroconversion_date))**alpha))
    error_term = random.gauss(0, eb*(Z**ea) + ec*(Z) + ed)

    if options.restrictBMV:
        return max(0, Z + error_term) 
    else:
        return Z + error_term
    #}}}

def ff3_cohort_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        param_set, 
        version,
        biol_id,
        fun = True): #{{{

    return bioObj.load_set(biol_id, version, param_set)[1] #}}}

def ff3_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id, 
        fun = True): #{{{

    alpha_mu = ff_cohort_pars['alpha_mu']
    alpha_sd = ff_cohort_pars['alpha_sd']
    alpha_trunc = ff_cohort_pars['alpha_trunc']
    beta_mu = ff_cohort_pars['beta_mu']
    beta_sd = ff_cohort_pars['beta_sd']
    beta_trunc = ff_cohort_pars['beta_trunc']
    height_mu = ff_cohort_pars['height_mu']
    height_sd = ff_cohort_pars['height_sd']
    height_trunc = ff_cohort_pars['height_trunc']

    def alpha_gen():
        x = -1000000
        while x < alpha_trunc:
            x = random.gauss(alpha_mu, alpha_sd)
        return x

    def beta_gen():
        x = -1000000
        while x < beta_trunc:
            x = random.gauss(beta_mu, beta_sd)
        return x

    def height_gen():
        x = -1000000
        while x < height_trunc:
            return random.gauss(height_mu, height_sd)

    results = odict()
    results['alpha'] = alpha_gen()
    results['beta'] = beta_gen()
    results['height'] = height_gen()

    return results #}}}

def ff3_exact():
    pass

biology3 = Biology(id=2, version="v1.0", table=bioObj, parameters = [#{{{
    AnnotatedTextSetting(name="biol2_par_alpha_mu", title="Alpha: Mean", default = 1.6143, db_var_name = "alpha_mu"),
    AnnotatedTextSetting(name="biol2_par_alpha_sd", title="Alpha: Standard Deviation", default = 0.2, db_var_name = "alpha_sd"),
    AnnotatedTextSetting(name="biol2_par_alpha_trunc", title="Alpha: Lower Limit", default = 0.1, db_var_name = "alpha_trunc"),
    AnnotatedTextSetting(name="biol2_par_beta_mu", title="Beta: Mean", default = 5, db_var_name = "beta_mu"),
    AnnotatedTextSetting(name="biol2_par_beta_sd", title="Beta: Standard Deviation", default = .5, db_var_name = "beta_sd"),
    AnnotatedTextSetting(name="biol2_par_beta_trunc", title="Beta: Lower Limit", default = 0.0001, db_var_name = "beta_trunc"),
    AnnotatedTextSetting(name="biol2_par_height_mu", title="Height: Mean", default = 1.95990, db_var_name = "height_mu"),
    AnnotatedTextSetting(name="biol2_par_height_sd", title="Height: Standard Deviation", default = 0.34214, db_var_name = "height_sd"),
    AnnotatedTextSetting(name="biol2_par_height_trunc", title="Height: Lower Limit", default = 0.5, db_var_name = "height_trunc"),
    AnnotatedTextSetting(name="biol2_par_ea", title="ea", default = 0, db_var_name = "ea"),
    AnnotatedTextSetting(name="biol2_par_eb", title="eb", default = 0, db_var_name = "eb"),
    AnnotatedTextSetting(name="biol2_par_ec", title="ec", default = 0, db_var_name = "ec"),
    AnnotatedTextSetting(name="biol2_par_ed", title="ed", default = 0.15, db_var_name = "ed"),
    ]) 

biology3['bmf_fun'] = ff3_bmf_fun
biology3['cohort_pars_gen'] = ff3_cohort_pars_gen
biology3['sub_pars_gen'] = ff3_sub_pars_gen
biology3['exact'] = ff3_exact
biology3['version'] = 'v1.0'
biology3['biol_id'] = 2 #}}}

# }}}

# Specify biology4 {{{ CD4 count for gbot
def ff4_bmf_fun(ff_cohort_pars, 
        ff_sub_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        visits): #{{{

    def epsilon(cd4):
        """
        COV as desribed in Hughes et al 1994 - ask Phillip
        """
        return 0.93 - 0.110 * math.log( max(1,cd4) )

    t = visits[-1][1] # Date of visit
    icd4 = ff_sub_pars['icd4']
    cd4dpy = ff_sub_pars['cd4dpy']
    arv_start = ff_sub_pars['arv_start']

    if t > arv_start:
        cd4dpy = -0.5 * cd4dpy

    Z = icd4 - (t/365.0) * cd4dpy
    cov = epsilon(tbmv)
    error_term = random.gauss(0, cov * tbmv)
    
    if options.restrictBMV:
        return max(0, Z + error_term) 
    else:
        return Z + error_term
    #}}}

def ff4_cohort_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        param_set, 
        version,
        biol_id,
        fun = True): #{{{

    return bioObj.load_set(biol_id, version, param_set)[1]
#    param_sets = odict()
#
#    combo1 = odict()
#    combo1['biol_param_set'] = 1
#    combo1['mu_picd4'] = 250
#    combo1['sd_picd4'] = 30
#    combo1['mu_cd4dpy'] = 60
#    combo1['sd_cd4dpy'] = 6
#    combo1['arv_start_l'] = 500
#    combo1['arv_start_u'] = 1000
#    combo1['arv_prob'] = 0.5
#    param_sets[1] = combo1
#
#    combo2 = odict()
#    combo2['biol_param_set'] = 2
#    combo2['mu_picd4'] = 250
#    combo2['sd_picd4'] = 30
#    combo2['mu_cd4dpy'] = 60
#    combo2['sd_cd4dpy'] = 6
#    combo2['arv_start_l'] = 500
#    combo2['arv_start_u'] = 1000
#    combo2['arv_prob'] = 0.1
#    param_sets[2] = combo2
#
#    return param_sets[param_set] #}}}

def ff4_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id, 

        fun = True): #{{{

    mu_picd4 = ff_cohort_pars['mu_picd4']
    sd_picd4 = ff_cohort_pars['sd_picd4']
    mu_cd4dpy = ff_cohort_pars['mu_cd4dpy']
    sd_cd4dpy = ff_cohort_pars['sd_cd4dpy']
    arv_start_l = ff_cohort_pars['arv_start_l'] 
    arv_start_u = ff_cohort_pars['arv_start_u'] 
    arv_prob  = ff_cohort_pars['arv_prob'] 


    def icd4_gen():
        return random.gauss(mu_picd4, sd_picd4)

    def cd4dpy_gen():
        return random.gauss(mu_cd4dpy, sd_cd4dpy)

    def arv_start_gen():
        if random.uniform(0,1) < arv_prob:
            return random.uniform(arv_start_l, arv_start_u)
        else:
            return 100000000

    results = odict()
    results['icd4'] = icd4_gen()
    results['cd4dpy'] = cd4dpy_gen()
    results['arv_start'] = arv_start_gen()

    return results #}}}

biology4 = Biology(id=3, version="v0.1", table=bioObj, parameters = [ #{{{
        AnnotatedTextSetting(name="biol3_par_mu_picd4", title="Initial CD4 mu", default = 250, db_var_name = "mu_picd4"),
        AnnotatedTextSetting(name="biol3_par_sd_picd4", title="Initial CD4 sigma", default = 30, db_var_name = "sd_picd4"),
        AnnotatedTextSetting(name="biol3_par_mu_cd4dpy", title="CD4 Decline PY mu", default = 60, db_var_name = "mu_cd4dpy"),
        AnnotatedTextSetting(name="biol3_par_sd_cd4dpy", title="CD4 Decline PY sigma", default = 6, db_var_name = "sd_cd4dpy"),
        AnnotatedTextSetting(name="biol3_par_arv_start_l", title="ARV initiation interval lower bound", default = 500, db_var_name = "arv_start_l"),
        AnnotatedTextSetting(name="biol3_par_arv_start_u", title="ARV initiation interval upper bound", default = 1000, db_var_name = "arv_start_u"),
        AnnotatedTextSetting(name="biol3_par_arv_prob", title = "ARV initiation probability", default = 0.5, db_var_name = "arv_prob")
        ]) 
biology4['bmf_fun'] = ff4_bmf_fun
biology4['cohort_pars_gen'] = ff4_cohort_pars_gen
biology4['sub_pars_gen'] = ff4_sub_pars_gen
biology4['version'] = 'v0.1'
biology4['biol_id'] = 3 #}}}

# }}}

# Specify ff5 {{{ specify distribution of crossing times and realize it with straight lines
def ff5_bmf_fun(ff_cohort_pars, 
        ff_sub_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        visits): #{{{

    t = visits[-1][1]
    sigma = ff_cohort_pars['sigma']
    slope = ff_sub_pars['slope']
    intercept = ff_sub_pars['intercept']
    
    error_term = random.gauss(0, sigma)
    Z = intercept + slope * t
    
    if options.restrictBMV:
        return max(0, Z + error_term) 
    else:
        return Z + error_term
    #}}}

def ff5_cohort_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        param_set, 
        version,
        biol_id,
        fun = True): #   {{{

    return bioObj.load_set(biol_id, version, param_set)[1] #}}}

def ff5_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id, 
        fun = True): #{{{

    dist = ff_cohort_pars['dist']
    dist_params = ff_cohort_pars['dist_params']
    sc = prot_sub_pars['seroconversion_date'] 

    if dist == "normal":
        eval_string = "random.gauss({dp})".format(dp = dist_params)

    slope = 1/float(eval(eval_string))
    intercept = -slope*sc

    results = odict()
    results['slope'] = slope
    results['intercept'] = intercept

    return results #}}}

biology5 = Biology(id = 4, 
        version = "v0.1", 
        table = bioObj, 
        parameters = [#{{{
            AnnotatedTextSetting(name="biol4_par_dist", 
                title="Name of Distribution", 
                default = "normal", 
                db_var_name = "dist"),
            AnnotatedTextSetting(name="biol4_par_dist_params", 
                title="Distribution parameters", 
                default = "180, 36", 
                db_var_name = "dist_params"),
            AnnotatedTextSetting(name="biol4_par_sigma", 
                title="Sigma", 
                default = 0.05, 
                db_var_name = "sigma"),
            AnnotatedTextSetting(name="biol4_par_assay_threshold", 
                title="Assay Threshold", 
                default = 1, 
                db_var_name = "assay_threshold")
            ]) 

biology5['bmf_fun'] = ff5_bmf_fun
biology5['cohort_pars_gen'] = ff5_cohort_pars_gen
biology5['sub_pars_gen'] = ff5_sub_pars_gen
biology5['version'] = "v0.1"
biology5['biol_id'] = 4 #}}}

## }}}

# Specify ff6 {{{ 4 parameter logistic
def ff6_bmf_fun(ff_cohort_pars, 
        ff_sub_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        visits): #{{{

    t = visits[-1][1]
    alpha = ff_sub_pars['alpha']
    beta = ff_sub_pars['beta']
    gamma = ff_sub_pars['gamma']
    delta = ff_sub_pars['delta']
    e0 = ff_cohort_pars['e0']
    e1 = ff_cohort_pars['e1']
    e2 = ff_cohort_pars['e2']
    e3 = ff_cohort_pars['e3']
    seroconversion_date = prot_sub_pars['seroconversion_date']

    Z = ((alpha - delta)/(1+(((t-seroconversion_date)/gamma)**beta)))+delta
    error_term = random.gauss(0, e2*(Z**e3) + e1*(Z) + e0)

    if options.restrictBMV:
        return max(0, Z + error_term) 
    else:
        return Z + error_term
    #}}}

def ff6_cohort_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        param_set, 
        version,
        biol_id,
        fun = True): #{{{

    return bioObj.load_set(biol_id, version, param_set)[1] #}}}

def ff6_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        prot_sub_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id, 
        fun = True): #{{{

    alpha_mu = ff_cohort_pars['alpha_mu']
    alpha_sd = ff_cohort_pars['alpha_sd']
    alpha_beta_sd = ff_cohort_pars['alpha_beta_sd']
    alpha_gamma_sd = ff_cohort_pars['alpha_gamma_sd']
    alpha_delta_sd = ff_cohort_pars['alpha_delta_sd']
    beta_mu = ff_cohort_pars['beta_mu']
    beta_sd = ff_cohort_pars['beta_sd']
    beta_gamma_sd = ff_cohort_pars['beta_gamma_sd']
    beta_delta_sd = ff_cohort_pars['beta_delta_sd']
    gamma_mu = ff_cohort_pars['gamma_mu']
    gamma_sd = ff_cohort_pars['gamma_sd']
    gamma_delta_sd = ff_cohort_pars['gamma_delta_sd']
    delta_mu = ff_cohort_pars['delta_mu']
    delta_sd = ff_cohort_pars['delta_sd']

    mymeans = [alpha_mu, beta_mu, gamma_mu, delta_mu]
    
    mycov = [[alpha_sd**2,       alpha_beta_sd,   alpha_gamma_sd,  alpha_delta_sd], 
            [alpha_beta_sd,   beta_sd**2,         beta_gamma_sd,   beta_delta_sd], 
            [alpha_gamma_sd,  beta_gamma_sd,   gamma_sd**2,        gamma_delta_sd],
            [alpha_delta_sd,  beta_delta_sd,   gamma_delta_sd,  delta_sd**2]]

    def accRejSam():
        done = False
        while not done:
            x = np.random.multivariate_normal(mymeans, mycov,1)
            if (x[0][1]>0) and (x[0][2]>0) and (x[0][3]>0) and (x[0][3]>x[0][0]):
                done = True
        return x

    alpha, beta, gamma, delta = list(accRejSam()[0])

    results = odict()
    results['alpha'] = alpha
    results['beta'] = beta
    results['gamma'] = gamma
    results['delta'] = delta

    return results #}}}

def ff6_exact(biol_id, version, param_set, method = "quad", threshold = 0, bigT = 365):

    def evalBMFnoNoise(t, beta, gamma, delta, seroconversion_date = 0):
        print beta, gamma, delta
        try:
            x = ((0 - delta)/(1+(((max(0,t-seroconversion_date))/gamma)**beta)))+delta
        except:
            print t, beta, gamma, delta
            raise Exception("your in trouble")
        return x

    def compNoise(Z, e0, e1, e2, e3):
        return e2*(Z**e3) + e1*(Z) + e0

    def comp_overthreshprob(t, beta, gamma, delta, seroconversion_date, e0, e1, e2, e3, threshold):
        Z = max(0, evalBMFnoNoise(t, beta, gamma, delta))
        return scipy.stats.norm.cdf((threshold - Z)/compNoise(Z, e0, e1, e2, e3))

    def comp_mdrpercurve(beta, gamma, delta, seroconversion_date, e0, e1, e2, e3, bigT, threshold):
        return scipy.integrate.quad(comp_overthreshprob, 0, bigT, args = (beta, gamma, delta, seroconversion_date, e0, e1, e2, e3, threshold))

    def dmvnorm(b, mean_vec, cov_mat):
        k = b.shape[0]
        part1 = numpy.exp(-0.5*k*numpy.log(2*numpy.pi))
        part2 = numpy.power(numpy.linalg.det(cov_mat),-0.5)
        dev = b-mean_vec
        part3 = numpy.exp(-0.5*numpy.dot(numpy.dot(dev.transpose(),numpy.linalg.inv(cov_mat)),dev))
        return part1*part2*part3

    def comp_mdrpercurve_pdf(delta, gamma, beta, seroconversion_date, e0, e1, e2, e3, bigT, threshold, mean_vec, cov_mat):
        b = np.array([beta, gamma, delta])
        print b
        #print mean_vec
        #print cov_mat
        #print "density"
        #print dmvnorm(b, mean_vec, cov_mat)
        return scipy.integrate.quad(comp_overthreshprob, 0, bigT, args = (beta, gamma, delta, seroconversion_date, e0, e1, e2, e3, threshold))[0] * dmvnorm(b, mean_vec, cov_mat)

    print method
    print threshold
    print bigT

    seroconversion_date = 0
    ff_cohort_pars = bioObj.load_set(biol_id, version, param_set)[1]
    e0 = ff_cohort_pars['e0']
    e1 = ff_cohort_pars['e1']
    e2 = ff_cohort_pars['e2']
    e3 = ff_cohort_pars['e3']
    the_score = 0
    old_score = 0
    start_time = time.time()
    if method == "montecarlo": #{{{
        step = 5000
        i = 0
        while True:
            i += 1
            ff_sub_pars = ff6_sub_pars_gen(0, 0, 0, ff_cohort_pars, 0, 0)
    
            alpha = ff_sub_pars['alpha']
            beta = ff_sub_pars['beta']
            gamma = ff_sub_pars['gamma']
            delta = ff_sub_pars['delta']
            print beta, gamma, delta
            sol = comp_mdrpercurve(beta, gamma, delta, seroconversion_date, e0, e1, e2, e3, bigT, threshold)[0]
            
            the_score += sol
            if i % step == 0:
                print i, round(sol,3), round(old_score/float(i-step) - the_score/float(i),3) if i > step else 0, round(the_score/float(i),10), round(time.time() - start_time), the_score
                old_score = the_score + 0 #}}}
    elif method == "quad":
        def comp_allcurve():
            # par inits {{{
            alpha_mu = ff_cohort_pars['alpha_mu']
            alpha_sd = ff_cohort_pars['alpha_sd']
            alpha_beta_sd = ff_cohort_pars['alpha_beta_sd']
            alpha_gamma_sd = ff_cohort_pars['alpha_gamma_sd']
            alpha_delta_sd = ff_cohort_pars['alpha_delta_sd']
            beta_mu = ff_cohort_pars['beta_mu']
            beta_sd = ff_cohort_pars['beta_sd']
            beta_gamma_sd = ff_cohort_pars['beta_gamma_sd']
            beta_delta_sd = ff_cohort_pars['beta_delta_sd']
            gamma_mu = ff_cohort_pars['gamma_mu']
            gamma_sd = ff_cohort_pars['gamma_sd']
            gamma_delta_sd = ff_cohort_pars['gamma_delta_sd']
            delta_mu = ff_cohort_pars['delta_mu']
            delta_sd = ff_cohort_pars['delta_sd']

            mymeans = [alpha_mu, beta_mu, gamma_mu, delta_mu]
            
            mycov = [[alpha_sd**2,       alpha_beta_sd,   alpha_gamma_sd,  alpha_delta_sd], 
                    [alpha_beta_sd,   beta_sd**2,         beta_gamma_sd,   beta_delta_sd], 
                    [alpha_gamma_sd,  beta_gamma_sd,   gamma_sd**2,        gamma_delta_sd],
                    [alpha_delta_sd,  beta_delta_sd,   gamma_delta_sd,  delta_sd**2]]
            #}}}

            mean_vec = [beta_mu, gamma_mu, delta_mu]
            
            cov_mat = [[beta_sd**2,         beta_gamma_sd,   beta_delta_sd], 
                    [beta_gamma_sd,   gamma_sd**2,        gamma_delta_sd],
                    [beta_delta_sd,   gamma_delta_sd,  delta_sd**2]]

            print mean_vec
            print cov_mat
            
            print scipy.integrate.tplquad(comp_mdrpercurve_pdf, 
                    max(0.1,beta_mu - 5*beta_sd), beta_mu + 5*beta_sd, 
                    lambda x: max(0.1,gamma_mu - 5*gamma_sd), lambda x: gamma_mu +5*gamma_sd,
                    lambda x,y : max(0.1,delta_mu - 5*delta_sd), lambda x,y : delta_mu + 5*delta_sd,
                    args = (seroconversion_date, e0, e1, e2, e3, bigT, threshold, mean_vec, cov_mat))
            print time.time() - start_time
        comp_allcurve()

            #scipy.integrate.tplquad(bigmamma2, 0, 1, lambda x: 0, lambda x: 1, lambda x, y: 0, lambda x, y: 1, args = (5,))

#    print comp_mdrpercurve(alpha, beta, gamma, delta, 0, e0, e1, e2, e3, bigT = 365, threshold = 40)

biology6 = Biology(id=6, version="v1.1", table=bioObj, parameters = [#{{{
    AnnotatedTextSetting(name="biol6_par_alpha_mu", title="Alpha: Mean", default = 0, db_var_name = "alpha_mu"),
    AnnotatedTextSetting(name="biol6_par_alpha_sd", title="Alpha: Standard Deviation", default = 0, db_var_name = "alpha_sd"),
    AnnotatedTextSetting(name="biol6_par_beta_mu", title="Beta: Mean", default = 4, db_var_name = "beta_mu"),
    AnnotatedTextSetting(name="biol6_par_beta_sd", title="Beta: Standard Deviation", default = 1.4, db_var_name = "beta_sd"),
    AnnotatedTextSetting(name="biol6_par_gamma_mu", title="Gamma: Mean", default = 190, db_var_name = "gamma_mu"),
    AnnotatedTextSetting(name="biol6_par_gamma_sd", title="Gamma: Standard Deviation", default = 50, db_var_name = "gamma_sd"),
    AnnotatedTextSetting(name="biol6_par_delta_mu", title="Delta: Mean", default = 85, db_var_name = "delta_mu"),
    AnnotatedTextSetting(name="biol6_par_delta_sd", title="Delta: Standard Deviation", default = 7.5, db_var_name = "delta_sd"),
    AnnotatedTextSetting(name="biol6_par_al_be_cov", title="Cov: Alpha and Beta", default = 0, db_var_name = "alpha_beta_sd"),
    AnnotatedTextSetting(name="biol6_par_al_ga_cov", title="Cov: Alpha and Gamma", default = 0, db_var_name = "alpha_gamma_sd"),
    AnnotatedTextSetting(name="biol6_par_al_de_cov", title="Cov: Alpha and Delta", default = 0, db_var_name = "alpha_delta_sd"),
    AnnotatedTextSetting(name="biol6_par_be_ga_cov", title="Cov: Beta and Gamma", default = -28, db_var_name = "beta_gamma_sd"),
    AnnotatedTextSetting(name="biol6_par_be_de_cov", title="Cov: Beta and Delta", default = 3.15000, db_var_name = "beta_delta_sd"),
    AnnotatedTextSetting(name="biol6_par_ga_de_cov", title="Cov: Gamma and Delta", default = -45, db_var_name = "gamma_delta_sd"),
    AnnotatedTextSetting(name="biol6_par_e3", title="Error: Power", default = 0.5 , db_var_name = "e3"),
    AnnotatedTextSetting(name="biol6_par_e2", title="Error: Power. coef.", default = 0.3 , db_var_name = "e2"),
    AnnotatedTextSetting(name="biol6_par_e1", title="Error: Lin. coef.", default = 0 , db_var_name = "e1"),
    AnnotatedTextSetting(name="biol6_par_e0", title="Error: Intercept", default = 2 , db_var_name = "e0"),
    ]) 

biology6['bmf_fun'] = ff6_bmf_fun
biology6['cohort_pars_gen'] = ff6_cohort_pars_gen
biology6['sub_pars_gen'] = ff6_sub_pars_gen
biology6['exact'] = ff6_exact
biology6['version'] = 'v1.1'
biology6['biol_id'] = 6 #}}}

# }}}

#biologies = [biology1]
#}}}
biologiesD = {'biology0':biology1,
        'biology1':biology2,
        'biology2':biology3,
        'biology3':biology4,
        'biology4':biology5,
        'biology6':biology6}
#}}}

# Specify Protocols {{{
# Description {{{
# The protocol describes the more admin side of the trial / study
# Where the function form family describes how the biomarker evolves,
#   the protocol describes how many subjects are tested,
#   how often they are tested, the probabilities that they
#   will miss a visit and so forth.
# Protocols are handled similarly to the funcion forms.
# One specifies a number of functions that generate gaps between the visits,
#   the visits that are missed, when the subject is lost to follow up and so forth
# Those functions depend on certain parameters. One then specify the distributions
#   from which those parameters are sampled. Lastly one speciy the parameters for
#   those distributions. The combination of this information can then be used to
#   simulate a large number of protocols to allow one to explore the effect of
#   for example increasing the lost to follow-up on the accuracy of the estimates
# A protocol family is specified by writing for functions:
# px_prot_pars_gen
#   This function generates the protocol specific parameters.
#   At the time this document was written, only the number of cohorts to
#   be simulated for a protocol needs to be specified in this function
# px_cohort_pars_gen
#   This function specifies the parameters that will remain the same for 
#   an entire cohort. Typically these should be hand-coded values that are
#   picked so that the situation is realistic. If one so choses, one may
#   allow these values to be sampled from a distribution. Typical values that
#   are specified in this function is mean gap between followup visits or
#   the mean probability that a subject will miss a visit.
# px_sub_pars_gen
#   This function will sample values for a subject. For example,
#   This function will return the probability that a subject
#   misses a visit. Most of these values will be sampled by executing
#   a function (with proper parameters passed to the function)
#   that was written in the shared function pool
# px_vis_par_gen
#   This function will return that values that are specific to the next visit
#   Basically, it only returns:
#   1) Was this the subjects last visit?
#   2) If not, when is the next visit
# }}}

# The protocol shared function pool {{{

# Description {{{
# A bunch of the protocols use the same functions to, 
#   for example, generate the gaps between visits
# Those functions are stored here and then only assigned
#   new names in the functions that specify the protocols
# This is just a way to keep from copy and pasting to much
# }}}

# The subject level functions of the protocol shared function pool {{{
prot_sub = {}

def tmp(fpv):
    if fpv > 0:
        #return round(random.uniform(0, fpv - 1)) old rounding table
        return random.uniform(0, fpv - 1)
    else:
        print "fpv <= 0 ERROR!"
        print 1/0

prot_sub['sc_date_gen_1'] = tmp

def tmp(alpha, beta, ti1):
    return ti1*random.betavariate(alpha, beta)

prot_sub['sc_date_gen_2'] = tmp

def tmp(vmpnMu, vmpnSigma):
    """
    Generate the probability of missing a visit while
    Subject is HIV negative
    """
    return random.gauss(vmpnMu, vmpnSigma)

prot_sub['vmpn_gen_1'] = tmp

def tmp(vmppMu, vmppSigma):
    """
    Generate the probability of missing a visit while
    Subject is HIV positive
    """
    return random.gauss(vmppMu, vmppSigma)

prot_sub['vmpp_gen_1'] = tmp

def tmp(vgnMu, vgnSigma):
    """
    Generate the gap between subsequent visits while
    Subject is HIV negative
    """
    done = False
    while not done:
        xasdf = random.gauss(vgnMu, vgnSigma)
        #print vgnMu, vgnSigma
        #print xasdf, 1
        if xasdf >= 0.0:
            done = True
    return xasdf

prot_sub['vgn_gen_1'] = tmp

def tmp(fpv_cutoff, vmpn, vgnMu, vgnSigma):
    """
    First positive visit date of subject
    """
    vgn_gen = prot_sub['vgn_gen_1']
    fpv = vgn_gen(vgnMu, vgnSigma)
    while random.uniform(0,1) < vmpn: #Did he miss the visit?
        fpv += vgn_gen(vgnMu, vgnSigma)
#    while fpv > fpv_cutoff:
#        fpv = vgn_gen(vgnMu, vgnSigma)
#        while random.uniform(0,1) < vmpn: #Did he miss the visit?
#            fpv += vgn_gen(vgnMu, vgnSigma)
    return fpv

prot_sub['fpv_gen_1'] = tmp
# }}}

# The visit level functions of the protocol shared function pool {{{
prot_vis = {}

def tmp(vgpMu, vgpSigma):
    """
    Generate the gap between subsequent visits while
    Subject is HIV positive
    """
    done = False
    while not done:
        xsdfg = random.gauss(vgpMu, vgpSigma)
        #print xsdfg, 2
        if xsdfg >= 0.0:
            done = True
    return xsdfg

prot_vis['vgp_gen_1'] = tmp
#}}}

# }}}

# Specify protocol1 - Fixed Study Duration {{{ 
def p1_prot_pars_gen(): #{{{
    results = odict()
    results['prot_id'] = 0
    results['prot_version'] = 'v2.2'
    results['ncohorts'] = 3
    
    return results #}}}

def p1_cohort_pars_gen(prot_prot_pars, 
        param_set): #{{{
    """
    Function to generate the cohort level parameters of a protocol
    """

    return protObj.load_set(prot_prot_pars['prot_id'], prot_prot_pars['prot_version'], param_set)[1]
    #}}}

def p1_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id): #{{{

    TMin = prot_cohort_pars['TMin']
    TMax  = prot_cohort_pars['TMax']
    ltfup  = prot_cohort_pars['ltfup']
    alphasc  = prot_cohort_pars['alphasc']
    betasc  = prot_cohort_pars['betasc']
    vmpnAlpha  = prot_cohort_pars['vmpnAlpha']
    vmpnBeta  = prot_cohort_pars['vmpnBeta']
    vmppAlpha  = prot_cohort_pars['vmppAlpha']
    vmppBeta  = prot_cohort_pars['vmppBeta']
    vgnMu  = prot_cohort_pars['vgnMu']
    vgnSigma  = prot_cohort_pars['vgnSigma']
    vgpMu  = prot_cohort_pars['vgpMu']
    vgpSigma = prot_cohort_pars['vgpSigma']
    tau = prot_cohort_pars['tau']

    sc_date_gen = prot_sub['sc_date_gen_2']
    vgn_gen = prot_sub['vgn_gen_1']
    fpv_gen = prot_sub['fpv_gen_1']

    def T_gen():
        if random.random() < ltfup: # Is he lost to followup?
            return random.uniform(TMin, TMax)
        else:
            return TMax

    T = T_gen()
    try:
        vmpn = random.betavariate(vmpnAlpha, vmpnBeta)
    except ValueError:
        vmpn = 0

    try:
        vmpp = random.betavariate(vmppAlpha, vmppBeta)
    except ValueError:
        vmpp = 0

    fpv = fpv_gen(T, vmpn, vgnMu, vgnSigma)
    seroconversion_date = sc_date_gen(alphasc, betasc, fpv)

    results = odict()
    results['seroconversion_date'] = seroconversion_date
    results['T'] = T
    results['vmpp'] = vmpp
    results['vmpn'] = vmpn
    results['fpv'] = fpv

    return results #}}}

def p1_prot_visit_pars_gen(prot_prot_pars,  
        prot_cohort_pars, 
        ff_cohort_pars, 
        prot_sub_pars, 
        ff_sub_pars, 
        visits,
        cohort_id, 
        sub_id, 
        visit_id): #{{{
    """
    Generates the next visit given 17 tons of info
    """
    TMin = prot_cohort_pars['TMin']
    TMax = prot_cohort_pars['TMax']
    ltfup = prot_cohort_pars['ltfup']
    vgnMu = prot_cohort_pars['vgnMu']
    vgnSigma = prot_cohort_pars['vgnSigma']
    vgpMu = prot_cohort_pars['vgpMu']
    vgpSigma = prot_cohort_pars['vgpSigma']
    tau = prot_cohort_pars['tau']
    ninc = prot_cohort_pars['ninc']

    done = False
    fpv = prot_sub_pars['fpv']
    vmpp = prot_sub_pars['vmpp']
    T = prot_sub_pars['T']

    vgp_gen = prot_vis['vgp_gen_1']

    if visits == []:
        visit_date = fpv
        return [done, [visit_id, visit_date, None]]
    else:
        x = vgp_gen(vgpMu, vgpSigma)
        while (random.uniform(0,1) < vmpp*((1+tau)**min(ninc, len(visits)-2))) and (visits[-1][1]+x <= (T+fpv)):
            x += vgp_gen(vgpMu, vgpSigma)
        if visits[-1][1]+x > (T+fpv):
            visit_date = visits[-1][1]+x
            visit_status = "missed"
            return [True, [visit_id, visit_date, visit_status]]
        else:
            visit_date = visits[-1][1]+x
            return [done, [visit_id, visit_date, None]] #}}}

protocol1 = Protocol(id=0, version="v2.2", table=protObj, parameters = [#{{{
        AnnotatedTextSetting(db_var_name = "cohort_sizes", name = "prot0_cohort_sizes", title = "Cohort sizes", default = 50, db_var_type = "Integer"),
        AnnotatedTextSetting(db_var_name = "TMin",name = "prot0_Tmin", title = "T minimum", default=0),
        AnnotatedTextSetting(db_var_name = "TMax",name = "prot0_Tmax", title = "T maximum", default=730),
        AnnotatedTextSetting(db_var_name = "ltfup",name = "prot0_ltfup", title = "Lost to Follow-up Probability", default=0),
        AnnotatedTextSetting(db_var_name = "alphasc",name = "prot0_alphasc", title = "Seroconversion Interval shape - Alpha", default=1),
        AnnotatedTextSetting(db_var_name = "betasc",name = "prot0_betasc", title = "Seroconversion Interval shape - Beta", default=1),
        AnnotatedTextSetting(db_var_name = "vmpnAlpha",name = "prot0_vmpnAlpha", title = "Visit Miss Probability Negative - Alpha", default=0),
        AnnotatedTextSetting(db_var_name = "vmpnBeta",name = "prot0_vmpnBeta", title = "Visit Miss Probability Negative - Beta", default=1),
        AnnotatedTextSetting(db_var_name = "vmppAlpha",name = "prot0_vmppAlpha", title = "Visit Miss Probability Positive - Alpha", default=0),
        AnnotatedTextSetting(db_var_name = "vmppBeta",name = "prot0_vmppBeta", title = "Visit Miss Probability - Beta", default=1),
        AnnotatedTextSetting(db_var_name = "vgnMu",name = "prot0_vgnMu", title = "Visit Gap Negative - Mu", default=90),
        AnnotatedTextSetting(db_var_name = "vgnSigma",name = "prot0_vgnSigma", title = "Visit Gap Negative - Sigma", default=9),
        AnnotatedTextSetting(db_var_name = "vgpMu",name = "prot0_vgpMu", title = "Visit Gap Positive - Mu", default=30),
        AnnotatedTextSetting(db_var_name = "vgpSigma",name = "prot0_vgpSigma", title = "Visit Gap Positive - Sigma", default=3),
        AnnotatedTextSetting(db_var_name = "tau",name = "prot0_tau", title = "Miss Probability Increase per Visit - r", default=0),
        AnnotatedTextSetting(db_var_name = "ninc",name = "prot0_ninc", title = "Number of Time the Visit Miss Probability can Increase - Ninc", default=0),
        ])
protocol1['prot_pars_gen']       = p1_prot_pars_gen
protocol1['cohort_pars_gen']     = p1_cohort_pars_gen
protocol1['sub_pars_gen']        = p1_sub_pars_gen
protocol1['visit_pars_gen'] = p1_prot_visit_pars_gen
#}}}
#}}}

# Specify protocol2 - Fixed number of SCHEDULED visits {{{
def p2_prot_pars_gen(): #{{{
    results = odict()
    results['prot_id'] = 2
    results['prot_version'] = 'v0.2'
    results['ncohorts'] = 3
    
    return results #}}}

def p2_cohort_pars_gen(prot_prot_pars, 
        param_set): #{{{
    """
    Function to generate the cohort level parameters of a protocol
    """
    return protObj.load_set(prot_prot_pars['prot_id'], prot_prot_pars['prot_version'], param_set)[1] #}}}

def p2_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id): #{{{

    vgn_mu  = prot_cohort_pars['vgn_mu']
    vgn_sd  = prot_cohort_pars['vgn_sd']
    vmpn    = prot_cohort_pars['vmpn']
    n_bin   = prot_cohort_pars['n_bin']
    p_bin   = prot_cohort_pars['p_bin']

    sc_date_gen = prot_sub['sc_date_gen_1']

    def fpv_gen(vmpn, vgn_mu, vgn_sd):
        fpv_gap = random.gauss(vgn_mu, vgn_sd)
        while random.uniform(0,1) < vmpn:
            fpv_gap += random.gauss(vgn_mu, vgn_sd)
        return fpv_gap

    n_visits = numpy.random.binomial(n_bin, p_bin)
    fpv = fpv_gen(vmpn, vgn_mu, vgn_sd)
    seroconversion_date = sc_date_gen(fpv)

    results = odict()
    results['seroconversion_date'] = seroconversion_date
    results['fpv'] = fpv
    results['n_visits'] = n_visits

    return results #}}}

def p2_prot_visit_pars_gen(prot_prot_pars,  
        prot_cohort_pars, 
        ff_cohort_pars, 
        prot_sub_pars, 
        ff_sub_pars, 
        visits,
        cohort_id, 
        sub_id, 
        visit_id): #{{{
    """
    Generates the next visit given 17 tons of info
    """
    vgp_mu      = prot_cohort_pars['vgp_mu']
    vgp_sd      = prot_cohort_pars['vgp_sd']
    n_visits    = prot_sub_pars['n_visits']
    vmpp        = prot_cohort_pars['vmpp']

    done = False
    fpv = prot_sub_pars['fpv']

    vgp_gen = prot_vis['vgp_gen_1']

    if visits == []:
        #visit_date = round(fpv) old rounding code
        visit_date = fpv
        visit_status = None
    else:
        x = vgp_gen(vgp_mu, vgp_sd)
        #visit_date = round(visits[-1][1]+x) old rounding code
        visit_date = visits[-1][1]+x
        if random.uniform(0,1) < vmpp:
            visit_status = "missed"
        else:
            visit_status = None
        if len(visits) >= n_visits:
            done = True
    return [done, [visit_id, visit_date, visit_status]] #}}}

protocol2 = Protocol(id=2, version='v0.2', table = protObj, parameters = [
    AnnotatedTextSetting(name="prot1_cohort_sizes", title="Cohort sizes", default = 30, db_var_name = "cohort_sizes", db_var_type = "Integer"),
    AnnotatedTextSetting(name="prot1_n_bin", title = "n Binomial", default = 10, db_var_name = "n_bin"),
    AnnotatedTextSetting(name="prot1_p_bin",title = "p Binomial", default = 0.9, db_var_name = "p_bin"),
    AnnotatedTextSetting(name="prot1_vmpn", title = "VMPN", default = 0.1, db_var_name = "vmpn"),
    AnnotatedTextSetting(name="prot1_vmpp", title = "VMPP", default = 0.1, db_var_name = "vmpp"),
    AnnotatedTextSetting(name="prot1_vgn_mu", title = "VGN mu", default = 90, db_var_name = "vgn_mu"),
    AnnotatedTextSetting(name="prot1_vgn_sd", title = "VGN sd", default = 5, db_var_name = "vgn_sd"),
    AnnotatedTextSetting(name="prot1_vgp_mu", title = "VGP mu", default = 30, db_var_name = "vgp_mu"),
    AnnotatedTextSetting(name="prot1_vgp_sd", title = "VGP sd", default = 3, db_var_name =  "vgp_sd")
    ])
protocol2['prot_pars_gen']       = p2_prot_pars_gen
protocol2['cohort_pars_gen']     = p2_cohort_pars_gen
protocol2['sub_pars_gen']        = p2_sub_pars_gen
protocol2['visit_pars_gen'] = p2_prot_visit_pars_gen
   #}}}

# Specify protocol3 - Exit when biomarker > threshold {{{
def p3_prot_pars_gen(): #{{{
    results = odict()
    results['prot_id'] = 1
    results['prot_version'] = 'v1.2'
    results['ncohorts'] = 3
    
    return results #}}}

def p3_cohort_pars_gen(prot_prot_pars, 
        param_set): #{{{
    """
    Function to generate the cohort level parameters of a protocol
    """
    return protObj.load_set(prot_prot_pars['prot_id'], prot_prot_pars['prot_version'], param_set)[1]
    #}}}

def p3_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id): #{{{

    vmpnMu  = prot_cohort_pars['vmpnMu']
    vmpnSigma  = prot_cohort_pars['vmpnSigma']
    vmppMu  = prot_cohort_pars['vmppMu']
    vmppSigma  = prot_cohort_pars['vmppSigma']
    vgnMu  = prot_cohort_pars['vgnMu']
    vgnSigma  = prot_cohort_pars['vgnSigma']
    vgpMu  = prot_cohort_pars['vgpMu']
    vgpSigma = prot_cohort_pars['vgpSigma']
    threshold_mult = prot_cohort_pars['threshold_mult']

    fpv_cutoff = 300

    sc_date_gen = prot_sub['sc_date_gen_1']
    vmpn_gen = prot_sub['vmpn_gen_1']
    vmpp_gen = prot_sub['vmpp_gen_1']
    vgn_gen = prot_sub['vgn_gen_1']
    fpv_gen = prot_sub['fpv_gen_1']

    def last_visit_prob_gen():
        return random.gauss(0.15, 0.015)

    vmpn = vmpn_gen(vmpnMu, vmpnSigma)
    fpv = fpv_gen(fpv_cutoff, vmpn, vgnMu, vgnSigma)
    seroconversion_date = sc_date_gen(fpv)
    vmpp = vmpp_gen(vmppMu, vmppSigma)
    last_visit_prob = last_visit_prob_gen()

    results = odict()
    results['seroconversion_date'] = seroconversion_date
    results['last_visit_prob'] = last_visit_prob
    results['vmpn'] = vmpn
    results['vmpp'] = vmpp
    results['fpv'] = fpv

    return results #}}}

def p3_prot_visit_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        ff_cohort_pars, 
        prot_sub_pars, 
        ff_sub_pars, 
        visits,
        cohort_id, 
        sub_id, 
        visit_id): #{{{
    """
    Generates the next visit given 17 tons of info
    """
    vmpnMu  = prot_cohort_pars['vmpnMu']
    vmpnSigma  = prot_cohort_pars['vmpnSigma']
    vmppMu  = prot_cohort_pars['vmppMu']
    vmppSigma  = prot_cohort_pars['vmppSigma']
    vgnMu  = prot_cohort_pars['vgnMu']
    vgnSigma  = prot_cohort_pars['vgnSigma']
    vgpMu  = prot_cohort_pars['vgpMu']
    vgpSigma = prot_cohort_pars['vgpSigma']
    threshold_mult = prot_cohort_pars['threshold_mult']

    done = False
    fpv = prot_sub_pars['fpv']
    vmpp = prot_sub_pars['vmpp']
    last_visit_prob = prot_sub_pars['last_visit_prob']
    assay_threshold = ff_cohort_pars['assay_threshold']

    vgp_gen = prot_vis['vgp_gen_1']

    if visits == []:
        visit_date = fpv
        #visit_date = round(fpv) old rounding code
    else:
        if random.uniform(0,1) < last_visit_prob:
            done = True
        if visits[-1][2] > threshold_mult * assay_threshold:
            done = True
        x = vgp_gen(vgpMu, vgpSigma)
        while random.uniform(0,1) < vmpp:
            x += vgp_gen(vgpMu, vgpSigma)
        visit_date = visits[-1][1]+x
        #visit_date = round(visits[-1][1]+x) old rounding code
    return [done, [visit_id, visit_date, None]] #}}}

protocol3 = Protocol(id=1, version='v0.2', table=protObj, parameters=[])
#AnnotatedTextSetting(name = "prot2_cohort_sizes"], title = "Cohort sizes", defaul = 30, db_var_name = "cohort_sizes", db_var_type = "Integer")
#AnnotatedTextSetting(name = "prot2_vmpnMu", title = "VMPN mu" , default = 0.2, db_var_name = "vmpnMu")
#AnnotatedTextSetting(name = "prot2_vmpnSigma" = 0.02
#AnnotatedTextSetting(name = "prot2_vmppMu" = 0.1
#AnnotatedTextSetting(name = "prot2_vmppSigma" = 0.01
#AnnotatedTextSetting(name = "prot2_vgnMu" = 40
#AnnotatedTextSetting(name = "prot2_vgnSigma" = 4
#AnnotatedTextSetting(name = "prot2_vgpMu" = 20
#AnnotatedTextSetting(name = "prot2_vgpSigma" = 2
#AnnotatedTextSetting(name = "prot2_threshold_mult" = 1.25
protocol3['prot_pars_gen']       = p3_prot_pars_gen
protocol3['cohort_pars_gen']     = p3_cohort_pars_gen
protocol3['sub_pars_gen']        = p3_sub_pars_gen
protocol3['visit_pars_gen'] = p3_prot_visit_pars_gen #}}}

# Specify protocol4 - Fixed number of REALIZED visits {{{
def p4_prot_pars_gen(): #{{{
    results = odict()
    results['prot_id'] = 3
    results['prot_version'] = 'v0.1'
    results['ncohorts'] = 3
    
    return results #}}}

def p4_cohort_pars_gen(prot_prot_pars, 
        param_set): #{{{
    """
    Function to generate the cohort level parameters of a protocol
    """
    return protObj.load_set(prot_prot_pars['prot_id'], prot_prot_pars['prot_version'], param_set)[1] #}}}

def p4_sub_pars_gen(prot_prot_pars, 
        prot_cohort_pars, 
        ff_cohort_pars, 
        cohort_id, 
        sub_id): #{{{

    vgn_lb   = prot_cohort_pars['vgn_lb']
    vgn_ub   = prot_cohort_pars['vgn_ub']

    vgp_time_dep_lb    = prot_cohort_pars['vgp_time_dep_lb']
    vgp_time_dep_ub    = prot_cohort_pars['vgp_time_dep_ub']

    vgp_time_dep = random.uniform(vgp_time_dep_lb, vgp_time_dep_ub)

    sc_date_gen = prot_sub['sc_date_gen_1']

    def fpv_gen(vgn_lb, vgn_ub):
        fpv = random.uniform(vgn_lb, vgn_ub)
        return fpv

    fpv = fpv_gen(vgn_lb, vgn_ub)
    seroconversion_date = sc_date_gen(fpv)

    results = odict()
    results['seroconversion_date'] = seroconversion_date
    results['fpv'] = fpv
    results['vgp_time_dep'] = vgp_time_dep

    return results #}}}

def p4_prot_visit_pars_gen(prot_prot_pars,  
        prot_cohort_pars, 
        ff_cohort_pars, 
        prot_sub_pars, 
        ff_sub_pars, 
        visits,
        cohort_id, 
        sub_id, 
        visit_id): #{{{
    """
    Generates the next visit given 17 tons of info
    """
    vgp_mu          = prot_cohort_pars['vgp_mu']
    vgp_sd          = prot_cohort_pars['vgp_sd']
    vgp_time_dep    = prot_sub_pars['vgp_time_dep']
    n_visits        = prot_cohort_pars['n_visits']

    done = False
    fpv = prot_sub_pars['fpv']

#    vgp_gen = prot_vis['vgp_gen_1']
    def vgp_gen(mu, sd):
        """
        finds alpha and beta for gamma from mu and sigma
        sim from random gamma
        """
        beta = sd/float(mu)
        alpha = mu/float(beta)
        return random.gammavariate(alpha, beta)

    if visits == []:
        #visit_date = round(fpv) old rounding code
        visit_date = fpv
        visit_status = None
    else:
        visit_number = len(visits)
        visit_time_dep_factor = vgp_time_dep**visit_number
        x = vgp_gen(vgp_mu * visit_time_dep_factor, vgp_sd * visit_time_dep_factor)
        visit_date = visits[-1][1]+x
        #visit_date = round(visits[-1][1]+x) # old rounding code
        visit_status = None
        if len(visits) >= (n_visits-1):
            done = True
    return [done, [visit_id, visit_date, visit_status]] #}}}

protocol4 = Protocol(id=3, version='v0.1', table = protObj, parameters = [
    AnnotatedTextSetting(name="prot4_cohort_sizes", title="Cohort sizes", default = 10, db_var_name = "cohort_sizes", db_var_type = "Integer"),
    AnnotatedTextSetting(name="prot4_vgp_mu", title = "Visit Gap while Positive: Mu", default = 50, db_var_name = "vgp_mu"),
    AnnotatedTextSetting(name="prot4_vgp_sd",title = "Visit Gap while Positive: Sigma", default = 10, db_var_name = "vgp_sd"),
    AnnotatedTextSetting(name="prot4_vgp_time_dep_lb", title = "Time Dependance of visit gaps: Lower Bound", default = 1.2, db_var_name = "vgp_time_dep_lb"),
    AnnotatedTextSetting(name="prot4_vgp_time_dep_ub", title = "Time Dependance of visit gaps: Upper Bound", default = 1.8, db_var_name = "vgp_time_dep_ub"),
    AnnotatedTextSetting(name="prot4_n_visits", title = "Number of Realized Visits", default = 5, db_var_name = "n_visits"),
    AnnotatedTextSetting(name="prot4_vgn_lb", title = "Visit Gap while Negative: Lower Bound", default = 30, db_var_name =  "vgn_lb"),
    AnnotatedTextSetting(name="prot4_vgn_ub", title = "Visit Gap while Positive: Upper Bound", default = 80, db_var_name =  "vgn_ub")
    ])
protocol4['prot_pars_gen']       = p4_prot_pars_gen
protocol4['cohort_pars_gen']     = p4_cohort_pars_gen
protocol4['sub_pars_gen']        = p4_sub_pars_gen
protocol4['visit_pars_gen'] = p4_prot_visit_pars_gen
   #}}}

protocolsD = {'prot0':protocol1,
#        'prot1':protocol3,
        'prot3':protocol4,
        'prot2':protocol2}
# }}}

# Functions to save data in the database {{{

def createTable(name, header): # {{{
    """
    This function will create a table that can contain the data
    """
    qstring = """
    CREATE TABLE %s
    (
    """ %(name,)
    for i in header:
        if i.count("id") == 1:
            colType = "int"
        elif (i.count("version") == 1) or (i == "prot_param_set") or (i == "bio_param_set"):
            colType = "varchar(50)"
        else:
            colType = "double"
        qstring += "%s %s,\n" %(i, colType,)
    qstring = qstring[:-2]
    qstring += ")Engine = InnoDB"
    return qstring #}}}

def dontQuoteNulls(x):
    if x == "null":
        return 'NULL'
    else:
        return "'"+str(x)+"'"

def insertData(name, data): #{{{
    qstring = """
    INSERT INTO %s VALUES (%s)
    """%(name, ",".join([dontQuoteNulls(bam) for bam in data]))
    return qstring #}}}

def roundVisitDates(name, roundit = True): # {{{
    if roundit:
        return """update %s set visit_date = round(visit_date,0)"""%name
    else:
        return """update %s set visit_date = visit_date"""%name
    #}}}

#}}}

def simulateCohorts(biologies, biolParSets, protocols, protParSets, ncohorts_input = -1): # {{{
    """
    Simulates all the data and saves it in MySQL
    """
    runid = "run_" + datetime.datetime.now().strftime("%Y_%m_%d_%H_%M") + "_" + "".join([str(random.randint(0,9)) for i in range(4)])
    print "runid:", runid
    for biology in biologies:
        bmf = biology['bmf_fun']
        ff_cohort_pars_gen = biology['cohort_pars_gen']
        ff_sub_pars_gen = biology['sub_pars_gen']
        biol_id = biology['biol_id']
        biol_version = biology['version']
        print "!"*50
        print "biol", biol_id
        print "!"*50
    
        for protocol in protocols:
            prot_prot_pars_gen = protocol['prot_pars_gen']
            prot_cohort_pars_gen = protocol['cohort_pars_gen']
            prot_sub_pars_gen = protocol['sub_pars_gen']
            prot_visit_pars_gen = protocol['visit_pars_gen']
    
            prot_prot_pars = prot_prot_pars_gen()
            prot_id = prot_prot_pars['prot_id']
            prot_version = prot_prot_pars['prot_version']
    
            print "-"*50
            print "prot", prot_id
            print "-"*50
    
            cohortTab = False
            subTab = False
            visitTab = False

            if ncohorts_input != -1:
                prot_prot_pars['ncohorts'] = ncohorts_input
                print "ncohorts overwritten with user input"
    
            protTabName = runid+"_prottab_" + str(biol_id) + "_" + str(prot_id)
            qstring = createTable(protTabName, 
                    ["biol_id"] + prot_prot_pars.keys() + ['biol_version'])
            con.execute(qstring)
            protTab = True
            qstring = insertData(name = protTabName,
                    data = [biol_id] + prot_prot_pars.values() + [biol_version])
            con.execute(qstring)
            db.commit()

            cohort_id = -1
            for prot_par_set in protParSets[prot_id]:
                for biol_par_set in biolParSets[biol_id]:
                    for cohort_counter in range(prot_prot_pars['ncohorts']):
                        cohort_id += 1
                        print "cohort", cohort_id, "  -><-  ",
                        
                        prot_cohort_pars = prot_cohort_pars_gen(prot_prot_pars, prot_par_set)
            
                        ff_cohort_pars = ff_cohort_pars_gen(prot_prot_pars, prot_cohort_pars, biol_par_set, biol_version, biol_id)
            
                        if not cohortTab:
                            cohortTabName = runid+"_cohorttab_" + str(biol_id) + "_" +str(prot_id)
                            qstring = createTable(cohortTabName, 
                                    ["biol_id", "prot_id", "cohort_id"] + prot_cohort_pars.keys() + ff_cohort_pars.keys())
                            con.execute(qstring)
                            cohortTab = True
                        #print [biol_id, prot_id, cohort_id] + prot_cohort_pars.values() + ff_cohort_pars.values()
                        qstring = insertData(name = cohortTabName,
                                data = [biol_id, prot_id, cohort_id] + prot_cohort_pars.values() + ff_cohort_pars.values())
                        con.execute(qstring)
                        db.commit()
#                        f = open('/tmp/sim2nd.log.tmp','a')
#                        f.writelines(str(prot_cohort_pars)+'\n')
#                        f.close()
                        
                        for sub_id in range(prot_cohort_pars['cohort_sizes']):
                            
                            prot_sub_pars = prot_sub_pars_gen(prot_prot_pars, prot_cohort_pars, ff_cohort_pars, cohort_id, sub_id)
            
                            ff_sub_pars = ff_sub_pars_gen(prot_prot_pars, prot_cohort_pars, prot_sub_pars, ff_cohort_pars, cohort_id, sub_id)
                            
                            if not subTab:
                                subTabName = runid+"_subtab_" + str(biol_id) + "_" + str(prot_id)
                                qstring = createTable(subTabName, 
                                        ["biol_id", "prot_id", "cohort_id", "sub_id"] + prot_sub_pars.keys() + ff_sub_pars.keys())
                                con.execute(qstring)
                                subTab = True
                            qstring = insertData(name = subTabName,
                                    data = [biol_id, prot_id, cohort_id, sub_id] + prot_sub_pars.values() + ff_sub_pars.values())
                            con.execute(qstring)
             
                            done = False
                            visit_id = -1
                            visits = []
                            while not done:
                                if not visitTab:
                                    visitTabName = runid+"_visittab_" + str(biol_id) + "_" + str(prot_id)
                                    qstring = createTable(visitTabName, 
                                            ["biol_id", "prot_id", "cohort_id", "sub_id", "visit_id", "visit_date", "bmv"])
                                    con.execute(qstring)
                                    visitTab = True
                                visit_id += 1
                                tmp = prot_visit_pars_gen(prot_prot_pars, 
                                        prot_cohort_pars, 
                                        ff_cohort_pars, 
                                        prot_sub_pars, 
                                        ff_sub_pars,
                                        visits, 
                                        cohort_id, 
                                        sub_id, 
                                        visit_id)
#                                f = open("/tmp/mysql_sim_2nd.log","a")
#                                f.writelines(str(tmp)+'\n')
#                                f.close()
                                done = tmp[0]
                                tmp_visits = visits + [tmp[1]]
                                if tmp[1][2] == "missed":
                                    pass
                                else:
                                    bmv = bmf(ff_cohort_pars, 
                                            ff_sub_pars, 
                                            prot_cohort_pars, 
                                            prot_sub_pars, 
                                            tmp_visits)
                                    tmp[1][2] = bmv
                                    visits.append(tmp[1])
                
                                    qstring = insertData(name = visitTabName,
                                            data = [biol_id, prot_id, cohort_id, sub_id]+tmp[1])
                                    #f = open("/tmp/mysql_sim_2nd.log","a")
                                    #f.writelines(qstring+'\n')
                                    #f.close()
                                    con.execute(qstring)
    # Only round visit dates if its not overwitten on the command line
    qstring = roundVisitDates(visitTabName, options.roundVisitDates)
    con.execute(qstring)
    db.commit()
    print "Success"
    return 0
    #}}}

def exactSolutions(biologies, biolParSets, method, threshold, bigT): #{{{
    """
    Computes the exact solutions for each biology in a list of biologies:
    """
    for biology in biologies:
        ff_cohort_pars_gen = biology['cohort_pars_gen']
        ff_exact = biology['exact']
        biol_id = biology['biol_id']
        biol_par_set = biolParSets[biol_id][0]
        ff_exact(biology['biol_id'], biology['version'], biol_par_set, method, threshold, bigT)

    # }}}

# prep protocols and biologies {{{
biologies = []
biolParamSets = {}
if options.reqBiologies != None:
    if options.reqBiologies.count(',') > 0:
        for i in options.reqBiologies.split(','):
            if i.count('.') == 0:
                biologies.append(biologiesD[i])
                biolParamSets[biologiesD[i]['biol_id']] = ['1']
            else:
                j = i.split(".")
                biologies.append(biologiesD[j[0]])
                biolParamSets[biologiesD[j[0]]['biol_id']] = []
                for k in j[1:]:
                    biolParamSets[biologiesD[j[0]]['biol_id']].append(k)
    else:
        i = options.reqBiologies
        if i.count('.') == 0:
            biologies.append(biologiesD[i])
            biolParamSets[biologiesD[i]['biol_id']] = ['1']
        else:
            j = i.split(".")
            biologies.append(biologiesD[j[0]])
            biolParamSets[biologiesD[j[0]]['biol_id']] = []
            for k in j[1:]:
                biolParamSets[biologiesD[j[0]]['biol_id']].append(k)

#print "biologies"
#print biologies
#print "biolParamSets"
#print biolParamSets

protocols = []
protParamSets = {}
if options.reqProtocols != None:
    if options.reqProtocols.count(',') > 0:
        for i in options.reqProtocols.split(','):
            if i.count('.') == 0:
                protocols.append(protocolsD[i])
                protParamSets[protocolsD[i]['prot_pars_gen']()['prot_id']] = ['1']
            else:
                j = i.split(".")
                protocols.append(protocolsD[j[0]])
                protParamSets[protocolsD[j[0]]['prot_pars_gen']()['prot_id']] = []
                for k in j[1:]:
                    protParamSets[protocolsD[j[0]]['prot_pars_gen']()['prot_id']].append(k)
    else:
        i = options.reqProtocols
        if i.count('.') == 0:
            protocols.append(protocolsD[i])
            protParamSets[protocolsD[i][0]()['prot_id']] = ['1']
        else:
            j = i.split(".")
            protocols.append(protocolsD[j[0]])
            protParamSets[protocolsD[j[0]]['prot_pars_gen']()['prot_id']] = []
            for k in j[1:]:
                protParamSets[protocolsD[j[0]]['prot_pars_gen']()['prot_id']].append(k)

#print "protocols"
#print protocols
#print "protParamSets"
#print protParamSets
# }}}

# Execute code {{{
if options.command == "exact":
    print('WARNING computation of exact solutions are buggy - rather use R')
    exactSolutions(biologies, biolParamSets, options.integration_method, options.threshold, options.bigT)
elif options.command == "sims":
    simulateCohorts(biologies, biolParamSets, protocols, protParamSets, int(options.ncohorts_input))
# }}}
