# -*- coding: utf-8 -*-
"""
Spyder Editor


author: Roberto Rosati
"""

from getpass import getpass
import requests
import json
import pandas as pd
import re
import numpy as np


# Helper functions
def notblank(info, secret = False):
    text = ''
    hidden = [input, getpass]
    while not text:
        text = hidden[secret](info.capitalize()+': ')
    return text

def printme(text):
    print(text)


def high_percent(number):
    return "{:.1%}".format(tofloat(number)/100)

def low_percent(number):
    return "{:.1%}".format(tofloat(number))

def flt(number):
    return "{:.1f}".format(tofloat(number))

def flt_M(number):
    return "{:.1f}M".format(tofloat(number)/1000000)
    

def integer(number):
    return "{:d}".format(number)

def higher(ref, val):
    return val >= ref

def lower(ref, val):
    return val <= ref

def tofloat(string):
    string = str(string)
    return float(string.replace('%', ''))


def pool(item):
    # Also works with the "old" notation with the pool at the end of the string
    pool = re.search('Pool=([\d]*?)(;|$)', item)
    if pool is not None:
        return int(pool.groups()[0])
    else:
        return False


# main class    
class Result(object):
    def __init__(self, server, auth, exp_id):
        printme("Analysing report ID: " + str(exp_id))
        printme("============================================================")
        self.data = {}
        if not server.startswith('http'):
            server = 'http://' + server
        while server.endswith('/'):
            server = server[:-1]
        self.server = server
        self.auth = auth
        self.exp_id = exp_id
        rel_url = '/rundb/api/v1/results/{}/'.format(self.exp_id)
        ret, ok = self.retrieve_json(rel_url)
        if ok:
            self.data['base_json'] = ret
            printme("Name:  " + self.data['base_json']['resultsName'])
            if self.data['base_json']['status'] != "Completed":
                printme('Status is not Completed: "{}"'.format(self.data['base_json']['status']))
                printme("Ending analysis.")
                return
            else:
                printme("Retrieving second links...")
                self.retrieve_datapages()
                printme("Checking plugin status...")
                ok = self.check_plugins()
                if ok:
                    printme('OK!')
                    ok = self.get_amplicon_tables()
                    if ok:
                        self.get_pool_data()
                        printme("Done!")
                        self.report()
                    else:
                        return
                else:
                    return
        else:
            return
        
    def retrieve_json(self, rel_url):
        abs_url = self.server+rel_url
        printme('...Retrieving: ' + abs_url)
        page = requests.get(abs_url, auth=auth, verify=False)
        if page.ok:
            out = json.loads(page.text)
        else:
            printme("(Error {}: {})".format(page.status_code, page.reason))
            out = None
        return out, page.ok
    
            
    
    def retrieve_datapages(self):
        for dataitem in ['analysismetrics',
                         'eas',
                         'experiment',
                         'libmetrics',
                         'pluginresults',
                         'qualitymetrics',
                         'resource_uri',
                         'tfmetrics']:
            this_object = self.data['base_json'][dataitem]
            if type(this_object) == str:
                ret, ok = self.retrieve_json(this_object)
                self.data[dataitem] = {'0': ret}
            elif type(this_object) == dict:
                self.data[dataitem] = {}
                for key, value in this_object:
                    self.data[dataitem][key] = self.retrieve_json(value)
            elif type(this_object) == list:
                self.data[dataitem] = {}
                for key, value in enumerate(this_object):
                    self.data[dataitem][key] = self.retrieve_json(value)
            else:
                printme("Unhandled object type: {}".format(type(this_object)))
                

    def check_plugins(self):
        self.plugin_IDs = {}
        for i in range(len(self.data['pluginresults'])):
            self.plugin_IDs[self.data['pluginresults'][i][0]['pluginName']] = i
            if self.data['pluginresults'][i][0]['state'] != 'Completed':
                printme("Plugin {} is not complete, returned status {}".format(
                            self.data['pluginresults'][i][0]['pluginName'],
                            self.data['pluginresults'][i][0]['state']))
                return False
        missing = set(['coverageAnalysis', 'variantCaller']).difference(set(self.plugin_IDs.keys()))
        if missing:
            printme("Results from these plugins were not found: {}".format(', '.join(missing)))
            return False
        return True
    
    def get_amplicon_tables(self):
        self.amplicon_tables = {}
        covan = self.data['pluginresults'][self.plugin_IDs['coverageAnalysis']][0]
        main_url = self.get_amplicons_main_url()
        if main_url:
            for barcode, jsondata in covan['store']['barcodes'].items():
                full_url = self.server + main_url + barcode + '/'+ jsondata['Alignments'] + '.amplicon.cov.xls'
                printme("Getting amplicon table for barcode: {}".format(barcode))
                data = requests.get(full_url, auth=self.auth, verify=False)
                data = [line.split('\t') for line in data.text.split('\n') if line != '']
                df = pd.DataFrame(data[1:], columns=data[0])
                df = df[['total_reads', 'attributes']]
                df['pool'] = df['attributes'].apply(pool)
                df = df[['total_reads', 'pool']]
                df['total_reads'] = df['total_reads'].astype(int)
                self.amplicon_tables[barcode] = df
            return True
        else:
            return False

    def get_amplicons_main_url(self):
        covan = self.data['pluginresults'][self.plugin_IDs['coverageAnalysis']][0]
        main_url = covan.get('URL', None)
        if not main_url:
            # fallback to 'path' for older versions of CoverageAnalysis
            main_url = covan.get('path', None)
            if main_url.startswith('/results/analysis'):
                main_url = main_url[17:] + '/'
            elif main_url is not None:
                printme("Unknown path format: " + main_url)
                main_url = None
        return main_url
                

    def get_pool_data(self):
        self.pools = {}
        for barcode in self.data['pluginresults'][self.plugin_IDs['coverageAnalysis']][0]['store']['barcodes']:
            self.pools[barcode] = []
            df = self.amplicon_tables[barcode]
            for pool in range(1, 13):
                self.pools[barcode].append(df[df['pool'] == pool]['total_reads'].mean())

    def get_stdev(self, barcode):
        dataset = self.pools[barcode]
        mean = tofloat(self.data['pluginresults'][self.plugin_IDs['coverageAnalysis']][0]['store']['barcodes'][barcode]['Average base coverage depth'])
        dataset = [100*datapoint/mean for datapoint in dataset]
        return float(np.std(dataset))
        
            
    def report(self):
        printme("============================================================")
        printme("Report ID:   " + str(self.data['base_json']['id']))
        printme("Report name: " + self.data['base_json']['resultsName'])
        printme("\nSamples:")
        samples = self.data['eas']['0']['barcodedSamples']
        for samplename, sampledata in samples.items():
            barcode = ','.join(list(sampledata['barcodeSampleInfo'].keys()))
            printme("  "+barcode+"  "+'"{}"'.format(samplename))
        printme("\nGlobal parameters:\n")
        # WARNING! For readability, the order of values and limits is inverted vs. the parameter table
        param_values = [("Loading", 80, self.data['analysismetrics'][0][0]['loading'], high_percent, higher),
                        ("Key signal", 70, self.data['libmetrics'][0][0]['aveKeyCounts'], integer, higher),
                        ("Mean raw accuracy", 98, self.data['libmetrics'][0][0]['raw_accuracy'], high_percent, higher)
                       ]
        self.parameter_table(param_values)
        printme("\nSample parameters:\n")
        samples = self.data['eas']['0']['barcodedSamples']
        for samplename, sampledata in samples.items():
            barcode = ','.join(list(sampledata['barcodeSampleInfo'].keys()))
            assert not ',' in barcode
            printme(barcode+"  "+'"{}"'.format(samplename))
            cov_an = self.data['pluginresults'][self.plugin_IDs['coverageAnalysis']][0]
            libmetrics = self.data['libmetrics'][0][0]
            # WARNING! For readability, the order of values and limits is inverted vs. the parameter table
            param_values = [("Number of mapped reads", 40000000, tofloat(cov_an['store']['barcodes'][barcode]['Number of mapped reads']), flt_M, higher),
                            ("Percent reads on target", 90, tofloat(cov_an['store']['barcodes'][barcode]['Percent reads on target']), high_percent, higher),
                            ("Average base coverage depth", 120, tofloat(cov_an['store']['barcodes'][barcode]['Average base coverage depth']), flt, higher),
                            ("Uniformity of base coverage", 90, tofloat(cov_an['store']['barcodes'][barcode]['Uniformity of base coverage']), high_percent, higher),
                            ("Percent Q20 bases", 0.8, tofloat(libmetrics['q20_mapped_bases'])/tofloat(libmetrics['q7_mapped_bases']), low_percent, higher),
                            ("Base coverage at 20x", 90, tofloat(cov_an['store']['barcodes'][barcode]['Target base coverage at 20x']), high_percent, higher),
                            ("Inter-pool standard dev.", 10, self.get_stdev(barcode), high_percent, lower)
                            ] 
            self.parameter_table(param_values)
            printme("")
        
        
            
            
                     
    def parameter_table(self, param_values):
        '''
        '''
        param_names = ['Parameter', 'Value', 'Limit', 'Outcome']
        col_widths = [30, # parameter name
                      6, # parameter value
                      6, # parameter limit
                      7 #fail/pass
                      ]
        dottedline = "-"*(sum(col_widths)+len(col_widths)+1)
        printme(dottedline)
        string = ' '
        for i, text in enumerate(param_names):
            string += text.ljust(col_widths[i]) + " "
        printme(string)
        printme(dottedline)
        for item in param_values:
            string = "|"
            string += item[0].ljust(col_widths[0]) + " "
            # Adjusting for the inverted order of values and limits
            for j, value in enumerate(item[2:0:-1]):
                string += item[3](value).rjust(col_widths[j+1]) + " "
            string += str(['FAIL', 'Pass'][int(item[4](item[1], item[2]))]).rjust(col_widths[3]) + '|'
            printme(string)
        printme(dottedline)
        

# might not be necessary
def search_result_json(server, auth, exp_id):
    if not server.startswith('http'):
        server = 'http://' + server

    results = server+'/rundb/api/v1/results/'
    
    my_exp_json_object = None
    # Check if experiment exists
    url = results
    while True:
        printme("Loading results page: " + url)
        page = requests.get(url, auth=auth, verify=False)
        page_json = json.loads(page.text)
        if 'objects' in page_json:
            for item in page_json['objects']:
                if item['id'] == exp_id:
                    my_exp_json_object = item
                    break
            if  my_exp_json_object:
                break
            elif page_json['meta']['next'] is not None:
                page = server + page_json['meta']['next']
                continue
            else:
                printme("Error: could not find report ID " + str(exp_id))
                break
        else:
            printme('Error: no objects retrieved from results API json')
            break
    return my_exp_json_object

    
    

if __name__ == '__main__':
    server = '172.16.178.1' #example
    user = 'ionadmin'
    auth = requests.auth.HTTPBasicAuth(user, notblank('Please input password for server "{}"'.format(server), secret=True))
    exp_id = int(notblank('Please input report ID number'))
    printme('')
    res = Result(server, auth, exp_id)
       
