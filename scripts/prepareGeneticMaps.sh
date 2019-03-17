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
GENMAPS='N.O.D.I.R'                ## -g
CHROMOSOME='N.O.C.H.R.O.M'         ## -c
OUTDIR='N.O.D.I.R'                 ## -o


################################################################################
# FUNCTIONS                                                                    #
################################################################################

# display_usage
# This function displays the usage of this program.
# No parameters
display_usage() {
  cat - <<EOF
  USAGE :
    ${NAME} [options] -g <in_genetic_maps_dir> -o <out_dir> -c <string>
      -g <inDirectory>  input directory where are stored genetic map files
      -c <string>       input chromosome
      -o <inDirectory>  directory where are stored output files
      -h                print help

  DESCRIPTION :
    ${NAME} removes header from genetic map file(s).

  EXAMPLE :
    ${NAME} -g genetic_map_HapMapII_GRCh37/ -o outdir/ -c 21
EOF

  return 0
}


################################################################################
# MAIN FUNCTION                                                                #
################################################################################
# main
# Removes header line from input genetic map files.
# Parameters : See 'getopts' part.
main() {

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
  while getopts :g:o:c: option
  do
    if [[ -z "${OPTARG}" ]]; then
      echo "[ERROR] Empty argument for option -${option}" >&2
      exit 1
    fi

    case "${option}" in
      g)
        GENMAPS="${OPTARG}"
        if [[ ! -d "${GENMAPS}" ]]; then
          echo "[ERROR] Input directory of genetic maps '${GENMAPS}' does not \
exist (option -g)." >&2
          exit 1
        fi
        ;; # -g <inDirectory>
      o)
        OUTDIR="${OPTARG}"
        if [[ ! -d "${OUTDIR}" ]]; then
          echo "[ERROR] Output directory '${OUTDIR}' does not exist \
(option -o). Please create it." >&2
          exit 1
        fi
        ;; # -o <inDirectory>
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

  readonly GENMAPS OUTDIR CHROMOSOME

  ### check mandatory parameters
  if [[ "${GENMAPS}" = 'N.O.D.I.R' ]]; then
    echo '[ERROR] Directory of genetic maps was not supplied (mandatory option \
-g)' >&2
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


  ### print used parameters
  cat - <<EOF
    ${NAME}
    Parameters as interpreted:
      -g ${GENMAPS}
      -c ${CHROMOSOME}
      -o ${OUTDIR}
EOF

  #fg_sar mark -l "prepGenMapsChr${CHROMOSOME}" || echo 'ignore fg_sar mark'

  #main process
  #because there is no genetic map for the Y chromosome, create empty files
  if [[ "${CHROMOSOME}" == "Y" ]]; then
    touch ${OUTDIR}/genetic_map_GRCh37_chr${CHROMOSOME}_wo_head.txt
  fi

  #because there is no genetic map for the Y chromosome
  if [[ "${CHROMOSOME}" != "Y" ]]; then

    if [[ "${CHROMOSOME}" != "X" ]]; then
      if [[ ! -f "${GENMAPS}/genetic_map_GRCh37_chr${CHROMOSOME}.txt" ]]; then
          echo "[ERROR] File '${GENMAPS}/genetic_map_GRCh37_chr${CHROMOSOME}.txt' does \
not exist. Program exit" >&2
          exit 1
      fi
      ## SORT by position (bp) - to get (lexically not numerically) sorted
      ## genetic maps as inputs of the future 'join' commands
      awk 'FNR>1 { $2=sprintf("%015d", $2); print $2,$4}' "${GENMAPS}/genetic_map_GRCh37_chr${CHROMOSOME}.txt" \
        | sort -k 1 \
        > "${OUTDIR}/genetic_map_GRCh37_chr${CHROMOSOME}_wo_head.txt"
      if [[ $? -ne 0 ]] \
        || [[ ! -f "${OUTDIR}/genetic_map_GRCh37_chr${CHROMOSOME}_wo_head.txt" ]] \
        || [[ ! -s "${OUTDIR}/genetic_map_GRCh37_chr${CHROMOSOME}_wo_head.txt" ]]; then
          echo "[ERROR] An error occurred when formating the genetic map file \
'${OUTDIR}/genetic_map_GRCh37_chr${CHROMOSOME}_wo_head.txt'. Program exit." >&2
          exit 1
      fi

    elif [[ "${CHROMOSOME}" == "X" ]]; then
      #there are 3 files for X
      if [[ ! -f "${GENMAPS}/genetic_map_GRCh37_chrX_par1.txt" ]]; then
        echo "[ERROR] File '${GENMAPS}/genetic_map_GRCh37_chrX_par1.txt' does not \
exist or is not a file. Program exit" >&2
        exit 1
      fi
      if [[ ! -f "${GENMAPS}/genetic_map_GRCh37_chrX.txt" ]]; then
        echo "[ERROR] File '${GENMAPS}/genetic_map_GRCh37_chrX.txt' does not \
exist or is not a file. Program exit" >&2
        exit 1
      fi
      if [[ ! -f "${GENMAPS}/genetic_map_GRCh37_chrX_par2.txt" ]]; then
        echo "[ERROR] File '${GENMAPS}/genetic_map_GRCh37_chrX_par2.txt' does not \
exist or is not a file. Program exit" >&2
        exit 1
      fi

      awk 'FNR>1 {print $2,$4}' "${GENMAPS}/genetic_map_GRCh37_chrX_par1.txt" \
        > "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt" \
        || {
          echo "[ERROR] An error occurred when formating the genetic map file \
'${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt'." >&2
          exit 1
      }
      awk 'FNR>1 {print $2,$4}' "${GENMAPS}/genetic_map_GRCh37_chrX.txt" \
        >> "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt" \
        || {
          echo "[ERROR] An error occurred when formating the genetic map file \
'${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt'." >&2
          exit 1
      }
      awk 'FNR>1 {print $2,$4}' "${GENMAPS}/genetic_map_GRCh37_chrX_par2.txt" \
        >> "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt"
      if [[ $? -ne 0 ]] \
        || [[ ! -f "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt" ]] \
        || [[ ! -s "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt" ]]; then
          echo "[ERROR] An error occurred when formating the genetic map file \
'${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt' or it is not a file." >&2
          exit 1
      fi

      ## SORT by position (bp) - to get (lexically not numerically) sorted
      ## genetic maps as inputs of the future 'join' commands
      awk '{
        $1=sprintf("%015d", $1);
        print $0
      }' "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt" \
        | sort -k 1 \
        > "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt.tmp" \
        || {
          echo "[ERROR] An error occurred when formating the genetic map file \
'${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt.tmp' or it is not a file." >&2
          exit 1
      }
      mv "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt.tmp" "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt"
      if [[ $? -ne 0 ]] \
        || [[ ! -f "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt" ]] \
        || [[ ! -s "${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt" ]]; then
          echo "[ERROR] An error occurred when formating the genetic map file \
'${OUTDIR}/genetic_map_GRCh37_chrX_wo_head.txt' or it is not a file." >&2
          exit 1
      fi
    fi
  fi

  return 0
}

main "$@"
