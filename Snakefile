#!/usr/bin/env python
# -*- coding: utf-8 -*-

################################################################################
# Copyright 2017 CEA CNRGH (Centre National de Recherche en Genomique Humaine) #
#                    <www.cnrgh.fr>                                            #
# Authors: Edith LE FLOCH (edith.le-floch@cea.fr)                              #
#          Elise LARSONNEUR (elise.larsonneur@cea.fr)                          #
#                                                                              #
# This software, LodSeq, is a computer program whose purpose is to perform     #
# genetic linkage analysis across families by computing lod-scores given a     #
# gVCF file and a related pedigree file.                                       #
#                                                                              #
# This software is governed by the CeCILL license under French law and         #
# abiding by the rules of distribution of free software.  You can  use,        #
# modify and/ or redistribute the software under the terms of the CeCILL       #
# license as circulated by CEA, CNRS and INRIA at the following URL            #
# "http://www.cecill.info".                                                    #
#                                                                              #
# As a counterpart to the access to the source code and  rights to copy,       #
# modify and redistribute granted by the license, users are provided only      #
# with a limited warranty  and the software's author,  the holder of the       #
# economic rights,  and the successive licensors  have only  limited           #
# liability.                                                                   #
#                                                                              #
# In this respect, the user's attention is drawn to the risks associated       #
# with loading,  using,  modifying and/or developing or reproducing the        #
# software by the user in light of its specific status of free software,       #
# that may mean  that it is complicated to manipulate,  and  that  also        #
# therefore means  that it is reserved for developers  and  experienced        #
# professionals having in-depth computer knowledge. Users are therefore        #
# encouraged to load and test the software's suitability as regards their      #
# requirements in conditions enabling the security of their systems and/or     #
# data to be ensured and,  more generally, to use and operate it in the        #
# same conditions as regards security.                                         #
#                                                                              #
# The fact that you are presently reading this means that you have had         #
# knowledge of the CeCILL license and that you accept its terms.               #
#                                                                              #
################################################################################

import subprocess
import re


shell.executable("/bin/bash")
shell.prefix("set -eo pipefail; ")

#python functions
#================
#function to define optional parameters, a default value can be assigned to the parameter when it does not exist in the config file
def getOptionalParam(cfg, key, defaultValue=None):
     try:
         value = cfg[key]
         return value
     except:
         return defaultValue

#global variables - software versions
#====================================
PLINK_OUTPUT = subprocess.check_output(["plink", "--version"]).decode('utf-8')
PLINK_VERSION = re.search(r'PLINK v([0-9]\.[0-9]{1,2}[a-zA-Z]*)', PLINK_OUTPUT).group(1)

VCFTOOLS_OUTPUT = subprocess.check_output(["vcftools"]).decode('utf-8')
VCFTOOLS_VERSION = re.search(r'VCFtools \(v?([0-9]\.[0-9]\.[0-9]{1,2})\)', VCFTOOLS_OUTPUT).group(1)

MERLIN_VERSION = ''
try:
   subprocess.check_output(["merlin"])
except subprocess.CalledProcessError as e:
   MERLIN_OUTPUT = e.output.decode('utf-8')
   MERLIN_VERSION = re.search(r'MERLIN ([0-9]\.[0-9]\.[0-9]{1,2})', MERLIN_OUTPUT).group(1)


#global variables - inputs
#=========================
chromosomes = config["chromosomes"]            ## if 'chromosomes' is a list in config file


#global variables - final outputs
#================================

#output files of rules runSinglePointMerlin and runMultiPointMerlin
OUTSGLDOMWOHTXT = expand('{out_dir}/runSinglePointMerlin/{ch}/results_singlepoint_chr{ch}_dominant.woheader.txt',
                out_dir=config["out_dir"],
                ch=chromosomes )

OUTSGLRECWOHTXT = expand('{out_dir}/runSinglePointMerlin/{ch}/results_singlepoint_chr{ch}_recessive.woheader.txt',
                out_dir=config["out_dir"],
                ch=chromosomes )

OUTMULTIDOMWOHTXT = expand('{out_dir}/runMultiPointMerlin/{ch}/results_multipoint_chr{ch}_dominant.woheader.txt',
                out_dir=config["out_dir"],
                ch=chromosomes )

OUTMULTIRECWOHTXT = expand('{out_dir}/runMultiPointMerlin/{ch}/results_multipoint_chr{ch}_recessive.woheader.txt',
                out_dir=config["out_dir"],
                ch=chromosomes )

OUTSGL = expand('{out_dir}/runSinglePointMerlin/{ch}/runSinglePointMerlin.out',
                out_dir=config["out_dir"],
                ch=chromosomes )

OUTMULTI = expand('{out_dir}/runMultiPointMerlin/{ch}/runMultiPointMerlin.out',
                out_dir=config["out_dir"],
                ch=chromosomes )

#output files of rule mergeResults
OUTMERGEDOMSGL = '{out_dir}/mergeResults/results_singlepoint_merged_dominant.txt'.format(out_dir=config["out_dir"])
OUTMERGERECSGL = '{out_dir}/mergeResults/results_singlepoint_merged_recessive.txt'.format(out_dir=config["out_dir"])
OUTMERGEDOMMULTI = '{out_dir}/mergeResults/results_multipoint_merged_dominant.txt'.format(out_dir=config["out_dir"])
OUTMERGERECMULTI = '{out_dir}/mergeResults/results_multipoint_merged_recessive.txt'.format(out_dir=config["out_dir"])


#localrules
#==========
# rules that are directly run by the main snakemake process
localrules: all


#onsuccess/onerror messages
##=========================
#Both handlers expect arbitrary python code (similar to the run keyword) and will be executed at the end of the workflow.
#As the names suggest, the first is executed if everything worked nicely, the second in case an error occurred.
#Handlers must be specified above rule definition
onstart:
    print("LodSeq version 1.0.0")
onsuccess:
    print("LodSeq finished with success.")
onerror:
    print("LodSeq failed.")
    #shell("echo 'See log file '{log}") #this log tmp file is removed...
    #shell("cat {log}")


#main rule contains target files
#main rule
#=========
rule all:
    input:
          [ OUTMERGEDOMSGL, OUTMERGERECSGL, OUTMERGEDOMMULTI, OUTMERGERECMULTI ]


#rules
#=====

#prepare genetic maps
rule prepareGeneticMaps:
    input:
        genetic_maps = config["genetic_maps"],
        out_dir = config["out_dir"]
    params:
        chromosomes = "{chromosomes}",
        out_log = config["out_log"]
    output:
        out_gen_map = config["out_dir"] + '/prepareGeneticMaps/'  + "{chromosomes}" + '/' + 'genetic_map_GRCh37_chr'  + "{chromosomes}" + '_wo_head.txt'
    log:
        out = config["out_dir"] + '/prepareGeneticMaps/'  + "{chromosomes}" + '/prepareGeneticMaps.out',
        err = config["out_dir"] + '/prepareGeneticMaps/'  + "{chromosomes}" + '/prepareGeneticMaps.err',
        rc = config["out_dir"] + '/prepareGeneticMaps/'  + "{chromosomes}" + '/prepareGeneticMaps.rc'
    message:
        'Running prepareGeneticMaps on chromosome {chromosomes}.'
    shell:
        """
        (
        if [[ "{params.chromosomes}" != "Y" ]]; then  #because there is no genetic map for the Y chromosome
          ## '{workflow.basedir}' to get the absolute path of snakefile directory
          ## use the absolute path '{workflow.basedir}/scripts' or the relative path './scripts' ?
          {workflow.basedir}/scripts/prepareGeneticMaps.sh \
            -g {input.genetic_maps} \
            -o {input.out_dir}/prepareGeneticMaps/{params.chromosomes}/ \
            -c {params.chromosomes}
        elif [[ "{params.chromosomes}" == "Y" ]]; then
          touch {output.out_gen_map}
        fi
        ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
        """


#get .ped and .map files from .vcf and .tfam files
#   and check that all input files and directories exist
rule prepareVCF:
    input:
        vcf = config["vcf"],
        tfam = config["tfam"],
        dom_model = config["dom_model"],
        rec_model = config["rec_model"],
        genetic_maps = config["genetic_maps"],
        out_dir = config["out_dir"]
    params:
        lod_threshold = config["lod_threshold"],
        out_prefix = config["out_prefix"],
        threads = getOptionalParam(config, "threads", '1'),
        out_log = config["out_log"]
    output:
        out_map = config["out_dir"] + '/' + 'prepareVCF' + '/' + config["out_prefix"] + '_vcftools_filled.map',
        out_ped = config["out_dir"] + '/' + 'prepareVCF' + '/' + config["out_prefix"] + '_vcftools_filled.ped'
    log:
        out = config["out_dir"] + '/prepareVCF/prepareVCF.out',
        err = config["out_dir"] + '/prepareVCF/prepareVCF.err',
        rc = config["out_dir"] + '/prepareVCF/prepareVCF.rc'
    #only for version change tracking - not in order to display version of a tool
    version:
        VCFTOOLS_VERSION
    version:
        PLINK_VERSION
    shell:
        """
        (
        if [[ ! -d $(dirname {params.out_log}) ]]; then echo -e 'Directory of log file does not exist; Program exit' 1>&2; exit 1; fi;
        echo '' 1> {params.out_log};
        if [[ ! -d {input.genetic_maps} ]]; then echo -e "Directory {input.genetic_maps} does not exist. Program exit" 1>&2; exit 1; fi;
        if [[ ! -d {input.out_dir} ]]; then echo -e "Directory {input.out_dir} does not exist. Program exit" 1>&2; exit 1; fi;
        if [[ ! -e {input.dom_model} ]]; then echo -e "Dominant model {input.dom_model} does not exist. Program exit" 1>&2; exit 1; fi;
        if [[ ! -e {input.rec_model} ]]; then echo -e "Recessive model {input.rec_model} does not exist. Program exit" 1>&2; exit 1; fi;
        if ! [[ {params.lod_threshold} =~ ^[0-9]+\.?[0-9]*$ ]] || [[ $(echo {params.lod_threshold} '<='0 | bc -l) -eq 1 ]]; then \
          echo -e 'The lod-score threshold must be greater than 0.' 1>&2; exit 1; fi;
        {workflow.basedir}/scripts/prepareVCF.sh \
          -i {input.vcf} \
          -p {input.tfam} \
          -o {input.out_dir}/prepareVCF/ \
          -s {params.out_prefix}_vcftools_filled \
          -t {params.threads}
        ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc}) ;
        """


#split files (.ped .map) by chromosome
#  then copy input files used by several rules to independent directories linked to each of this rule
#  to avoid the use of same input files by 2 tasks simultaneously
rule splitByChrom:
    input:
        out_dir = config["out_dir"],
        in_map = config["out_dir"] + '/' + 'prepareVCF' + '/' + config["out_prefix"] + '_vcftools_filled.map',
        in_ped = config["out_dir"] + '/' + 'prepareVCF' + '/' + config["out_prefix"] + '_vcftools_filled.ped'
    params:
        chromosomes = "{chromosomes}",
        out_prefix = config["out_prefix"],
        out_log = config["out_log"],
        threads = getOptionalParam(config, "threads", '1')
    threads:
        2
    output:
        out_chr_map=config["out_dir"] + '/splitByChrom/' + "{chromosomes}" + '/' + config["out_prefix"] + '_vcftools_filled_chr' + "{chromosomes}" + '.map',
        out_chr_ped=config["out_dir"] + '/splitByChrom/' + "{chromosomes}" + '/' + config["out_prefix"] + '_vcftools_filled_chr' + "{chromosomes}" + '.ped',
    log:
        out = config["out_dir"] + '/splitByChrom/'  + "{chromosomes}" + '/splitByChrom.out',
        err = config["out_dir"] + '/splitByChrom/'  + "{chromosomes}" + '/splitByChrom.err',
        rc = config["out_dir"] + '/splitByChrom/'  + "{chromosomes}" + '/splitByChrom.rc'
    shell:
        """
        (
        #without renaming outputs of plink at the end of plink, ie outputs of plink are expected as output files of snakemake process => successfull in mode local but NOT in mode cluster
        /bin/bash -ce '{workflow.basedir}/scripts/splitByChrom.sh \
          -m {input.in_map} \
          -p {input.in_ped} \
          -c {params.chromosomes} \
          -o {input.out_dir}/splitByChrom/{params.chromosomes}/ \
          -s {params.out_prefix}_vcftools_filled_chr \
          -t {params.threads};
        ';
        ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
        """


#prepare single point merlin analysis inputs
rule prepareSinglePointFiles:
    input:
        out_dir = config["out_dir"],
        in_chr_map = config["out_dir"] + '/splitByChrom/' + "{chromosomes}" + '/' + config["out_prefix"] + '_vcftools_filled_chr' + "{chromosomes}" + '.map',
    params:
        chromosomes = "{chromosomes}",
        out_prefix = config["out_prefix"],
        out_log = config["out_log"]
    output:
        out_chr_dat = config["out_dir"] + '/prepareSinglePointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_sgl_chr' + "{chromosomes}" + '.dat',
        out_chr_map = config["out_dir"] + '/prepareSinglePointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_sgl_chr' + "{chromosomes}" + '.map',
    log:
        out = config["out_dir"] + '/prepareSinglePointFiles/'  + "{chromosomes}" + '/prepareSinglePointFiles.out',
        err = config["out_dir"] + '/prepareSinglePointFiles/'  + "{chromosomes}" + '/prepareSinglePointFiles.err',
        rc = config["out_dir"] + '/prepareSinglePointFiles/'  + "{chromosomes}" + '/prepareSinglePointFiles.rc'
    shell:
        """
        (
        #outputs .dat and .map files
        {workflow.basedir}/scripts/prepareSinglePointFiles.sh \
          -m {input.in_chr_map} \
          -c {params.chromosomes} \
          -o {input.out_dir}/prepareSinglePointFiles/{params.chromosomes}/ \
          -s {params.out_prefix}_sgl_chr
        #copy the ped file ouputted by splitByChrom
        ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
        """


#run merlin in a single point way
rule runSinglePointMerlin:
    input:
        out_dir = config["out_dir"],
        dom_model = config["dom_model"],
        rec_model = config["rec_model"],
        in_chr_dat = config["out_dir"] + '/prepareSinglePointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_sgl_chr' + "{chromosomes}" + '.dat',
        in_chr_map = config["out_dir"] + '/prepareSinglePointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_sgl_chr' + "{chromosomes}" + '.map',
        in_chr_ped = config["out_dir"] + '/splitByChrom/' + "{chromosomes}" + '/' + config["out_prefix"] + '_vcftools_filled_chr' + "{chromosomes}" + '.ped'
    params:
        chromosomes = "{chromosomes}",
        lod_threshold = config["lod_threshold"],
        out_prefix = config["out_prefix"],
        out_log = config["out_log"]
    output:
        out_dom_txt = config["out_dir"] + '/runSinglePointMerlin/' + "{chromosomes}" + '/results_singlepoint_chr' + "{chromosomes}" + '_dominant.txt',
        out_rec_txt = config["out_dir"] + '/runSinglePointMerlin/' + "{chromosomes}" + '/results_singlepoint_chr' + "{chromosomes}" + '_recessive.txt',
        out_dom_signif_txt = config["out_dir"] + '/runSinglePointMerlin/' + "{chromosomes}" + '/results_singlepoint_chr' + "{chromosomes}" + '_dominant_LODsignif.txt',
        out_rec_signif_txt = config["out_dir"] + '/runSinglePointMerlin/' + "{chromosomes}" + '/results_singlepoint_chr' + "{chromosomes}" + '_recessive_LODsignif.txt',
        out_dom_woheader_txt = config["out_dir"] + '/runSinglePointMerlin/' + "{chromosomes}" + '/results_singlepoint_chr' + "{chromosomes}" + '_dominant.woheader.txt',
        out_rec_woheader_txt = config["out_dir"] + '/runSinglePointMerlin/' + "{chromosomes}" + '/results_singlepoint_chr' + "{chromosomes}" + '_recessive.woheader.txt',
        out_log_out = config["out_dir"] + '/runSinglePointMerlin/'  + "{chromosomes}" + '/runSinglePointMerlin.out'
    log:
        out = config["out_dir"] + '/runSinglePointMerlin/'  + "{chromosomes}" + '/runSinglePointMerlin.out',
        err = config["out_dir"] + '/runSinglePointMerlin/'  + "{chromosomes}" + '/runSinglePointMerlin.err',
        rc = config["out_dir"] + '/runSinglePointMerlin/'  + "{chromosomes}" + '/runSinglePointMerlin.rc'
    version:
        MERLIN_VERSION
    shell:
        """
        (
        {workflow.basedir}/scripts/runSinglePointMerlin.sh \
          -D {input.dom_model} \
          -R {input.rec_model} \
          -m {input.in_chr_map} \
          -d {input.in_chr_dat} \
          -p {input.in_chr_ped} \
          -c {params.chromosomes} \
          -o {input.out_dir}/runSinglePointMerlin/{params.chromosomes}/ \
          -s results_singlepoint_chr \
          -l {params.lod_threshold}
        ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
        """


#prepare multipoint merlin analysis inputs
rule prepareMultiPointFiles:
    input:
        out_dir = config["out_dir"],
        in_chr_map = config["out_dir"] + '/splitByChrom/' + "{chromosomes}" + '/' + config["out_prefix"] + '_vcftools_filled_chr' + "{chromosomes}" + '.map',
        in_chr_ped = config["out_dir"] + '/splitByChrom/' + "{chromosomes}" + '/' + config["out_prefix"] + '_vcftools_filled_chr' + "{chromosomes}" + '.ped',
        in_chr_genmap = config["out_dir"] + '/prepareGeneticMaps/'  + "{chromosomes}" + '/' + 'genetic_map_GRCh37_chr'  + "{chromosomes}" + '_wo_head.txt'
    params:
        chromosomes = "{chromosomes}",
        out_prefix = config["out_prefix"],
        out_log = config["out_log"],
        threads = getOptionalParam(config, "threads", '1')
    output:
        out_chr_dat = config["out_dir"] + '/prepareMultiPointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_multi_chr' + "{chromosomes}" + '.dat',
        out_chr_map = config["out_dir"] + '/prepareMultiPointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_multi_chr' + "{chromosomes}" + '.map',
        out_chr_ped = config["out_dir"] + '/prepareMultiPointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_multi_chr' + "{chromosomes}" + '.ped'
    log:
        out = config["out_dir"] + '/prepareMultiPointFiles/'  + "{chromosomes}" + '/prepareMultiPointFiles.out',
        err = config["out_dir"] + '/prepareMultiPointFiles/'  + "{chromosomes}" + '/prepareMultiPointFiles.err',
        rc = config["out_dir"] + '/prepareMultiPointFiles/'  + "{chromosomes}" + '/prepareMultiPointFiles.rc'
    shell:
        """
        (
        GENMAPDIR=$(dirname {input.in_chr_genmap});
        {workflow.basedir}/scripts/prepareMultiPointFiles.sh \
          -m {input.in_chr_map} \
          -p {input.in_chr_ped} \
          -c {params.chromosomes} \
          -g ${{GENMAPDIR}} \
          -o {input.out_dir}/prepareMultiPointFiles/{params.chromosomes}/ \
          -s {params.out_prefix}_multi_chr \
          -t {params.threads}
        ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
        """


#run merlin in a multipoint way
rule runMultiPointMerlin:
    input:
        out_dir = config["out_dir"],
        dom_model = config["dom_model"],
        rec_model = config["rec_model"],
        in_chr_dat = config["out_dir"] + '/prepareMultiPointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_multi_chr' + "{chromosomes}" + '.dat',
        in_chr_map = config["out_dir"] + '/prepareMultiPointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_multi_chr' + "{chromosomes}" + '.map',
        in_chr_ped = config["out_dir"] + '/prepareMultiPointFiles/' + "{chromosomes}" + '/' + config["out_prefix"] + '_multi_chr' + "{chromosomes}" + '.ped'
    params:
        chromosomes = "{chromosomes}",
        lod_threshold = config["lod_threshold"],
        out_prefix = config["out_prefix"],
        out_log = config["out_log"]
    output:
        out_dom_txt = config["out_dir"] + '/runMultiPointMerlin/' + "{chromosomes}" + '/results_multipoint_chr' + "{chromosomes}" + '_dominant.txt',
        out_rec_txt = config["out_dir"] + '/runMultiPointMerlin/' + "{chromosomes}" + '/results_multipoint_chr' + "{chromosomes}" + '_recessive.txt',
        out_dom_signif_txt = config["out_dir"] + '/runMultiPointMerlin/' + "{chromosomes}" + '/results_multipoint_chr' + "{chromosomes}" + '_dominant_LODsignif.txt',
        out_rec_signif_txt = config["out_dir"] + '/runMultiPointMerlin/' + "{chromosomes}" + '/results_multipoint_chr' + "{chromosomes}" + '_recessive_LODsignif.txt',
        out_dom_woheader_txt = config["out_dir"] + '/runMultiPointMerlin/' + "{chromosomes}" + '/results_multipoint_chr' + "{chromosomes}" + '_dominant.woheader.txt',
        out_rec_woheader_txt = config["out_dir"] + '/runMultiPointMerlin/' + "{chromosomes}" + '/results_multipoint_chr' + "{chromosomes}" + '_recessive.woheader.txt',
        out_log_out = config["out_dir"] + '/runMultiPointMerlin/'  + "{chromosomes}" + '/runMultiPointMerlin.out'
    log:
        out = config["out_dir"] + '/runMultiPointMerlin/'  + "{chromosomes}" + '/runMultiPointMerlin.out',
        err = config["out_dir"] + '/runMultiPointMerlin/'  + "{chromosomes}" + '/runMultiPointMerlin.err',
        rc = config["out_dir"] + '/runMultiPointMerlin/'  + "{chromosomes}" + '/runMultiPointMerlin.rc'
    shell:
        """
        (
        if [[ "{params.chromosomes}" != "Y" ]]; then  #because there is no genetic map for the Y chromosome
          {workflow.basedir}/scripts/runMultiPointMerlin.sh \
            -D {input.dom_model} \
            -R {input.rec_model} \
            -m {input.in_chr_map} \
            -d {input.in_chr_dat} \
            -p {input.in_chr_ped} \
            -c {params.chromosomes} \
            -o {input.out_dir}/runMultiPointMerlin/{params.chromosomes}/ \
            -s results_multipoint_chr \
            -l {params.lod_threshold}
        elif [[ "{params.chromosomes}" == "Y" ]]; then
             echo "Ignoring merlin multipoint analysis of the chromosome {params.chromosomes}.";
             touch {output.out_dom_txt};
             touch {output.out_rec_txt};
             touch {output.out_dom_signif_txt};
             touch {output.out_rec_signif_txt};
             touch {output.out_dom_woheader_txt};
             touch {output.out_rec_woheader_txt};
        fi;
        ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
        """



# merge final lod-score files
rule mergeResults:
    input:
        out_dir = config["out_dir"],
        in_sgl_dom_txt = OUTSGLDOMWOHTXT,
        in_sgl_rec_txt = OUTSGLRECWOHTXT,
        in_multi_dom_txt = OUTMULTIDOMWOHTXT,
        in_multi_rec_txt = OUTMULTIRECWOHTXT,
        in_sgl_log = OUTSGL,
        in_multi_log = OUTMULTI
    params:
        out_log = config["out_log"]
    output:
        out_merge_sgl_dom_txt = OUTMERGEDOMSGL,
        out_merge_sgl_rec_txt = OUTMERGERECSGL,
        out_merge_multi_dom_txt = OUTMERGEDOMMULTI,
        out_merge_multi_rec_txt = OUTMERGERECMULTI
    log:
        out = config["out_dir"] + '/mergeResults/mergeResults.out',
        err = config["out_dir"] + '/mergeResults/mergeResults.err',
        rc = config["out_dir"] + '/mergeResults/mergeResults.rc'
    shell:
        """
        (
        cat {input.in_sgl_dom_txt} 1> {output.out_merge_sgl_dom_txt}
        cat {input.in_sgl_rec_txt} 1> {output.out_merge_sgl_rec_txt}
        cat {input.in_multi_dom_txt} 1> {output.out_merge_multi_dom_txt}
        cat {input.in_multi_rec_txt} 1> {output.out_merge_multi_rec_txt}
        echo -e '\nSINGLEPOINT ANALYSIS' >> {params.out_log}
        for f in {input.in_sgl_log}; do if grep -q 'max' $f ; then grep 'max' $f >> {params.out_log}; fi; done
        echo -e '\nMULTIPOINT ANALYSIS' >> {params.out_log}
        for f in {input.in_multi_log}; do if grep -q 'max' $f ; then grep 'max' $f >> {params.out_log}; fi; done
        ) 1> {log.out} 2> {log.err} && echo $? > {log.rc} || echo $? > {log.rc} ; exit $(cat {log.rc});
        """

