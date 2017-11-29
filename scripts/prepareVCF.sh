#!/bin/bash

################################################################################
# Copyright 2017 CEA CNRGH (Centre National de Recherche en Genomique Humaine) #
#                    <www.cnrgh.fr>                                            #
# Authors: Edith LE FLOCH (edith.le-floch@cea.fr)			       #
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


set -uo pipefail


# check bash version support of regex with capturing groups
bash -c '[[ "a/b" =~ ^([a-z]+)/([b-x]*)$ ]] && [ "${BASH_REMATCH[1]}" = a ] && true' \
  &> /dev/null \
    || {
    echo "$(basename "$0"): bash >= 3.1 required." >&2
    exit 1
}


################################################################################
# GLOBAL VARIABLES                                                             #
################################################################################

NAME="$(basename "$0")"
readonly NAME

## parameters initialization
THREADS=1                          ## -t // option '--threads' of plink is ignored by plink...
VCF='N.O.F.I.L.E'                  ## -i
TFAM='N.O.F.I.L.E'                 ## -p
OUTDIR='N.O.D.I.R'                 ## -o
OUTPREFIX='N.O.S.T.R.I.N.G'        ## -s


################################################################################
# FUNCTIONS                                                                    #
################################################################################

# display_usage
# This function displays the usage of this program.
# No parameters
function display_usage {
  cat - <<EOF
  USAGE :
    ${NAME} [options] -i <in_vcf> -p <in_tfam> -o <out_dir>
      -i <inFile>       input vcf variant file
      -p <inFile>       input tfam pedigree file
      -o <inDirectory>  directory where are stored output files
      -s <string>       prefix of output files
      -t <int>          number of threads used by plink steps (default : 1)
      -h                print help
  
  DESCRIPTION :
    ${NAME} outputs .ped and .map files from input .vcf and .tfam files.

  EXAMPLE :
    ${NAME} -i cohort.vcf -p pedigree.tfam -o outdir/ -s cohort_result -t 1
EOF
  return 0
}


################################################################################
# MAIN FUNCTION                                                                #
################################################################################
# main
# Outputs .ped and .map files from input .vcf and .tfam files.
# Parameters : See 'getopts' part.
function main {

  local nb_proc=1 
  local vcftools_tmp_log=""
  local plink_tmp_log=""
  local rc=0
  local grep_rc=1

  # check whether user had supplied -h or --help . If yes display usage
  if [[ "$#" -eq 0 ]]; then
    display_usage
    exit 0
  fi
  if [[ -n "$1" ]]; then
    if [[ "$1" = "-?" ]] || [[ "$1" = "-h" ]] || [[ "$1" = "--help" ]]; then
      display_usage
      exit 0
    fi
  fi

  # if less than three arguments supplied, display usage
  if [[  "$#" -lt 3 ]]; then
    echo '[ERROR] Missing argument(s)' >&2
    display_usage
    exit 1
  fi


 ## catch option values
  while getopts :i:p:o:s:t: option
  do
    if [[ -z "${OPTARG}" ]]; then 
      echo "[ERROR] Empty argument for option -${option}" >&2 
      exit 1 
    fi

    case "${option}" in
      i)
        VCF="${OPTARG}"
        if [[ ! -f "${VCF}" ]]; then 
          echo "[ERROR] Input VCF file '${VCF}' does not exist or is not a file \
(option -i)." >&2 
          exit 1  
        fi
        ;; # -i <inFile>
      p)
        TFAM="${OPTARG}";
        if [[ ! -f "${TFAM}" ]]; then 
          echo "[ERROR] Input TFAM file '${TFAM}' does not exist or is not a file \
(option -p)." >&2 
          exit 1  
        fi
        ;; # -p <inFile>
      o)
        OUTDIR="${OPTARG}";
        if [[ ! -d "${OUTDIR}" ]]; then 
          echo "[ERROR] Output directory '${OUTDIR}' does not exist (option -o). \
Please create it." >&2 
          exit 1  
        fi
        ;; # -o <inDirectory>
      s)
        OUTPREFIX="${OPTARG}"
        ;; # -s <string>
      t)
        THREADS="${OPTARG}"
        nb_proc="$(getconf _NPROCESSORS_ONLN)"  #unix and osx
        if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || [[ "${THREADS}" -lt 1 ]]; then 
          echo '[ERROR] the number of threads must be greater than 0 (option -t).' >&2 
          exit 1  
        fi
        if [[ "${THREADS}" -gt "${nb_proc}" ]]; then 
          echo "[ERROR] too much threads requested, use ${nb_proc} thread(s) instead" >&2
          exit 1 
        fi
        ;; # -t <number of threads>
      :)
        echo "[ERROR] option ${OPTARG} : missing argument" >&2
        exit 1
        ;;
      \?)
        echo "[ERROR] ${OPTARG} : invalid option" >&2
        exit 1
        ;;
    esac
  done

  readonly VCF TFAM OUTDIR THREADS

  ### checking input directories and files
  if [[ "${VCF}" = 'N.O.F.I.L.E' ]]; then 
    echo '[ERROR] Input VCF file was not supplied (mandatory option -i)' >&2 
    exit 1  
  fi
  if [[ "${TFAM}" = 'N.O.F.I.L.E' ]]; then 
    echo '[ERROR] Input TFAM file was not supplied (mandatory option -p)' >&2 
    exit 1  
  fi
  if [[ "${OUTDIR}" = 'N.O.D.I.R' ]]; then 
    echo '[ERROR] Output directory was not supplied (mandatory option -o)' >&2 
    exit 1  
  fi
 
  ### define default output file prefix 
  if [[ "${OUTPREFIX}" = 'N.O.S.T.R.I.N.G' ]]; then
    OUTPREFIX="$(basename "${VCF}" '.vcf')" 
  fi
  readonly OUTPREFIX


  ### print used parameters
  cat - <<EOF
  ${NAME} 
  Parameters as interpreted:
    -i ${VCF} \n  -p ${TFAM} \n  -o ${OUTDIR} \n  -s ${OUTPREFIX} \n  -t ${THREADS}
EOF


  #main process
  ##display available chromosomes
  echo ''
  echo '[INFO] chromosomes for which there are variants in the vcf file: '
  zcat "${VCF}" | \
    awk '{
      if ( $1 !~ /^#/ && !($1 in a)){ 
        a[$1] = $1 
      } 
    } END {
      for (i in a) {
        print a[i]}
    }' \
    | sort -u \
    || echo '[ERROR] Program failed to extract chromosomes from variant file.' >&2
  echo ''


  ##transform from (g)vcf to tped/tfam format
  vcftools_tmp_log="$(mktemp "${OUTDIR}/vcftools.XXX.log")"
  # echo "vcftools --gzvcf ${VCF} --plink-tped --out ${OUTDIR}/${OUTPREFIX}_vcftools" 
  vcftools --gzvcf "${VCF}" --plink-tped --out "${OUTDIR}/${OUTPREFIX}_vcftools" \
    &> "${vcftools_tmp_log}"
  ##=> generate a tfam file with values of 0
  rc=$?
  cat "${vcftools_tmp_log}" | tee >(grep -i 'error' >&2) | grep -vi 'error'
  grep -qi 'error' "${vcftools_tmp_log}"
  grep_rc=$?
  if [[ $rc -ne 0 ]] \
    || [[ ${grep_rc} -eq 0 ]] \
    || [[ ! -f "${OUTDIR}/${OUTPREFIX}_vcftools.tped" ]] \
    || [[ ! -s "${OUTDIR}/${OUTPREFIX}_vcftools.tped" ]]; then 
      if [[ -f "${vcftools_tmp_log}" ]]; then rm "${vcftools_tmp_log}"; fi
      echo "[ERROR] File '${OUTDIR}/${OUTPREFIX}_vcftools.tped' was not output by \
vcftools or is empty. See log files. Program exit" >&2 
      exit 1 
  fi
  if [[ -f "${vcftools_tmp_log}" ]]; then rm "${vcftools_tmp_log}"; fi


  ##transform from tped/tfam to ped/map
  plink_tmp_log="$(mktemp "${OUTDIR}/plink.XXX.log")"
  plink \
    --tped "${OUTDIR}/${OUTPREFIX}_vcftools.tped" \
    --tfam "${TFAM}" \
    --recode \
    --out "${OUTDIR}/${OUTPREFIX}" \
    --threads "${THREADS}" \
    &> "${plink_tmp_log}"
  rc=$?
  cat "${plink_tmp_log}" | tee >(sed -n '/^Error/,$p' >&2) | grep -v '^Error'
  if [[ -f "${plink_tmp_log}" ]]; then rm "${plink_tmp_log}"; fi
  if [[ $rc -ne 0 ]] \
    || [[ ! -f "${OUTDIR}/${OUTPREFIX}.ped" ]] \
    || [[ ! -s "${OUTDIR}/${OUTPREFIX}.ped" ]] \
    || [[ ! -f "${OUTDIR}/${OUTPREFIX}.map" ]] \
    || [[ ! -s "${OUTDIR}/${OUTPREFIX}.map" ]]; then
      echo "[ERROR] plink failed while generating the files \
'${OUTDIR}/${OUTPREFIX}.ped' and '${OUTDIR}/${OUTPREFIX}.map'. \
See log files for more details." >&2
      exit $rc
  fi

  return 0
}

main "$@"
