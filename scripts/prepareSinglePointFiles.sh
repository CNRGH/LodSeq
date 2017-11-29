#!/bin/bash

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
CHROMOSOME='N.O.C.H.R.O.M'        ## -c
OUTDIR='N.O.D.I.R'                ## -o
OUTPREFIX='N.O.S.T.R.I.N.G'       ## -s


################################################################################
# FUNCTIONS                                                                    #
################################################################################

# display_usage
# This function displays the usage of this program.
# No parameters
function display_usage {
  cat - <<EOF
  USAGE :
    ${NAME} [options] -m <in_map> -c <string> -o <out_dir> 
      -m <inFile>       input map file
      -c <string>       input chromosome
      -o <inDirectory>  directory where are stored output files
      -s <string>       prefix of output files
      -h                print help
  
  DESCRIPTION :
    ${NAME} formats .dat and .map files for a singlepoint merlin analysis from an
    input .map file.

  EXAMPLE :
    ${NAME} -m cohort.map -c 21 -o outdir/ -s cohort_sgl_chr
EOF

  return 0
}


################################################################################
# MAIN FUNCTION                                                                #
################################################################################
# main
# Formats .dat and .map files for a singlepoint merlin analysis from an input .map file.
# Parameters : See 'getopts' part.
function main {

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
  while getopts :m:c:o:s: option
  do
    if [[ -z "${OPTARG}" ]]; then 
      echo "[ERROR] Empty argument for option -${option}" >&2 
      exit 1 
    fi

    case "${option}" in
      m)
        MAP="${OPTARG}"
        if [[ ! -f "${MAP}" ]]; then 
          echo "[ERROR] Input MAP file '${MAP}' does not exist or is not a file \
(option -m)." >&2 
          exit 1  
        fi
        ;; # -m <inFile>
      c)
        CHROMOSOME="${OPTARG}"
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

  readonly MAP OUTDIR CHROMOSOME


  ### checking input directories and files
  if [[ "${MAP}" = 'N.O.F.I.L.E' ]]; then 
    echo '[ERROR] Input MAP file was not supplied (mandatory option -m)' >&2 
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
    OUTPREFIX="$(basename "${MAP}" '.map')_sgl" 
  fi
  readonly OUTPREFIX


  ### print used parameters
  cat - <<EOF
    ${NAME} 
    Parameters as interpreted:
      -m ${MAP}
      -c ${CHROMOSOME}
      -o ${OUTDIR}
      -s ${OUTPREFIX}

EOF


  #main process
  #make a .dat file for merlin with bi-allelic markers filtered in by vcftools
  echo 'A affection' > "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat" \
    && awk '{print "M",$2}' "${MAP}" >> "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat" 
  if [[ $? -ne 0 ]] \
    || [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat" ]] \
    || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat" ]]; then
      echo "[ERROR] File '${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat' was not output \
or is not a file or is empty." >&2
      exit 1
  fi

  #remove the forth column (position in bp) of the map file #for merlin
  awk 'BEGIN{OFS="\t"};{print $1,$2,$3}' "${MAP}" > "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map"
  if [[ $? -ne 0 ]] \
    || [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map" ]] \
    || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map" ]]; then
      echo "[ERROR] File '${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map' was not output \
or is not a file or is empty." >&2
      exit 1
  fi

  return 0  
}

main "$@"
