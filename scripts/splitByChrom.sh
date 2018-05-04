#!/usr/bin/env bash

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
MAP='N.O.F.I.L.E'                 ## -m
PED='N.O.F.I.L.E'                 ## -p
CHROMOSOME='N.O.C.H.R.O.M'        ## -c
OUTDIR='N.O.D.I.R'                ## -o
OUTPREFIX='N.O.S.T.R.I.N.G'       ## -s
THREADS=1                         ## -t // option '--threads' of plink is ignored by plink...


################################################################################
# FUNCTIONS                                                                    #
################################################################################

# display_usage
# This function displays the usage of this program.
# No parameters
function display_usage {
  cat - <<EOF
  USAGE :
    ${NAME} [options] -m <in_map> -p <in_ped> -c <string> -o <out_dir>
      -m <inFile>       input map file
      -p <inFile>       input ped pedigree file
      -c <string>       input chromosome
      -o <inDirectory>  directory where are stored output files
      -s <string>       prefix of output files
      -t <int>          number of threads used by plink steps (default : 1)
      -h                print help
  
  DESCRIPTION :
    ${NAME} outputs .map and .ped files of a specified chromosome from input .map 
    and .ped files.

  EXAMPLE :
    ${NAME} -m cohort.map -p pedigree.ped -c 21 -o outdir/ -s cohort_chr -t 1
EOF

  return 0
}


################################################################################
# MAIN FUNCTION                                                                #
################################################################################
# main
# Performs a genetic linkage analysis (lod-scores) between related individuals.
# Parameters : See 'getopts' part.
function main {

  local nb_proc=1
  local plink_tmp_log=""
  local rc=0

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

  # if less than four arguments supplied, display usage
  if [[  "$#" -lt 4 ]]; then
    echo '[ERROR] Missing argument(s)' >&2
    display_usage
    exit 1
  fi



  ## catch option values
  while getopts :m:p:c:o:s:t: option
  do
    if [[ -z "${OPTARG}" ]]; then 
      echo "[ERROR] Empty argument for option -${option}" >&2 
      exit 1 
    fi

    case "${option}" in
      m)
        MAP="${OPTARG}";
        if [[ ! -f "${MAP}" ]]; then 
          echo "[ERROR] Input MAP file '${MAP}' does not exist or is not a file\
 (option -m)." >&2 
          exit 1  
        fi
        ;; # -m <inFile>
      p)
        PED="${OPTARG}";
        if [[ ! -f "${PED}" ]]; then 
          echo "[ERROR] Input PED file '${PED}' does not exist or is not a file\
 (option -p)." >&2 
          exit 1 
        fi
        ;; # -p <inFile>
      c)
        CHROMOSOME="${OPTARG}";
        if [[ "${CHROMOSOME}" =~ ^[0-9]+$ ]]; then
          if [[ "${CHROMOSOME}" -lt 1 ]] || [[ "${CHROMOSOME}" -gt 22 ]]; then 
            echo "[ERROR] invalid chromosome '${CHROMOSOME}' (option -c)." >&2 
            exit 1 
          fi
        else
          if [[ "${CHROMOSOME}" != "X" ]] && [[ "${CHROMOSOME}" != "Y" ]]; then 
            echo "[ERROR] invalid chromosome '${CHROMOSOME}' (option -c)." >&2
            exit 1 
          fi
        fi
        ;; # -c <string>
      o)
        OUTDIR="${OPTARG}"
        if [[ ! -d "${OUTDIR}" ]]; then 
          echo "[ERROR] Output directory '${OUTDIR}' does not exist (option -o).\
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
          echo "[ERROR] too much threads requested, use ${nb_proc} thread(s) \
instead" >&2
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

  readonly MAP PED OUTDIR CHROMOSOME THREADS

  ### check mandatory parameters
  if [[ "${MAP}" = 'N.O.F.I.L.E' ]]; then 
    echo '[ERROR] Input MAP file was not supplied (mandatory option -m)' >&2
    exit 1 
  fi
  if [[ "${PED}" = 'N.O.F.I.L.E' ]]; then 
    echo '[ERROR] Input PED file was not supplied (mandatory option -p)' >&2 
    exit 1 
  fi
  if [[ "${OUTDIR}" = 'N.O.D.I.R' ]]; then 
    echo '[ERROR] Output directory was not supplied (mandatory option -o)' >&2 
    exit 1 
  fi
  if [[ "${CHROMOSOME}" = 'N.O.C.H.R.O.M' ]]; then
    echo '[ERROR] Chromosome was not supplied (mandatory option -c)' >&2 
    exit 1 
  fi

  ### define default output file prefix
  if [[ "${OUTPREFIX}" = 'N.O.S.T.R.I.N.G' ]]; then 
    OUTPREFIX="$(basename "${MAP}" '.map')_chr" 
  fi
  readonly OUTPREFIX


  ### print used parameters
  cat - <<EOF
    ${NAME} 
    Parameters as interpreted:
      -m ${MAP} 
      -p ${PED} 
      -c ${CHROMOSOME} 
      -o ${OUTDIR} 
      -s ${OUTPREFIX} 
      -t ${THREADS}

EOF


  #main process
  cp -p "${MAP}" "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map" || exit $?
  cp -p "${PED}" "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped" || exit $?
 
  plink_tmp_log="$(mktemp "${OUTDIR}/plink.XXX.log")" 
  plink \
    --file "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK" \
    --chr "${CHROMOSOME}" \
    --recode \
    --out "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}" \
    --threads "${THREADS}" \
    &> "${plink_tmp_log}"
  rc=$?
  cat "${plink_tmp_log}" | tee >(sed -n '/^Error/,$p' >&2) | grep -v '^Error'
  if [[ -f "${plink_tmp_log}" ]]; then rm "${plink_tmp_log}"; fi
  if [[ $rc -ne 0 ]] \
    || [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map" ]] \
    || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map" ]] \
    || [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.ped" ]] \
    || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.ped" ]]; then
      echo "[ERROR] plink failed while generating the files \
'${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map' \
'${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.ped'. \
See log files for more details." >&2
      #remove temp files
      if [[ -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map" ]]; then 
        rm "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map" 
      fi
      if [[ -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped" ]]; then 
        rm "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped"
      fi
      exit $rc
  fi

  #remove temp files
  if [[ -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map" ]]; then 
    rm "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map" 
  fi
  if [[ -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped" ]]; then 
    rm "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped"
  fi

  return 0
}

main "$@" 
