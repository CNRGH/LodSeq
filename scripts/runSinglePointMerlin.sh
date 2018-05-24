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
DOM_MODEL='N.O.F.I.L.E'           ## -D
REC_MODEL='N.O.F.I.L.E'           ## -R
MAP='N.O.F.I.L.E'                 ## -m
DAT='N.O.F.I.L.E'                 ## -d
PED='N.O.F.I.L.E'                 ## -p
CHROMOSOME='N.O.C.H.R.O.M'        ## -c
OUTDIR='N.O.D.I.R'                ## -o
OUTPREFIX='N.O.S.T.R.I.N.G'       ## -s
MINLODTH=2.0                      ## -l

END_WITH_ERRORS=0


################################################################################
# FUNCTIONS                                                                    #
################################################################################

# display_usage
# This function displays the usage of this program.
# No parameters
display_usage() {
  cat - <<EOF
  USAGE :
    ${NAME} [options] -D <in_dom_model> -R <in_rec_model> -m <in_map> -d <in_dat> \
-p <in_ped> -c <string> -o <out_dir>
      -D <inFile>       input merlin dominant model file
      -R <inFile>       input merlin recessive model file
      -m <inFile>       input map file
      -d <inFile>       input dat file
      -p <inFile>       input ped pedigree file
      -c <string>       input chromosome
      -o <inDirectory>  directory where are stored output files
      -s <string>       prefix of output files
      -l <float>        minimal lod-score threshold (default : 2.0)
      -h                print help

  DESCRIPTION :
    ${NAME} runs a singlepoint merlin analysis from input .map .dat .ped files.

  EXAMPLE :
    ${NAME} -D dominant.model -R recessive.model -m cohort_sgl_chr21.map \
-d cohort_sgl_chr21.dat -p cohort_sgl_chr21.ped -c 21 -o outdir/ \
-s results_sgl_chr -l 2.0
EOF

  return 0
}


################################################################################
# MAIN FUNCTION                                                                #
################################################################################
# main
# Runs a singlepoint merlin analysis from input .map .dat .ped files.
# Parameters : See 'getopts' part.
main() {

  local file=""
  local prefix=""
  local line_nb=
  local end_file=

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

  # if less than seven arguments supplied, display usage
  if [[  "$#" -lt 7 ]]; then
    echo '[ERROR] Missing argument(s)' >&2
    display_usage
    exit 1
  fi



  ## catch option values
  while getopts :D:R:m:d:p:c:o:s:l: option
  do
    if [[ -z "${OPTARG}" ]]; then
      echo "[ERROR] Empty argument for option -${option}" >&2
      exit 1
    fi

    case "${option}" in
      D)
        DOM_MODEL="${OPTARG}"
        if [[ ! -f "${DOM_MODEL}" ]]; then
          echo "[ERROR] Input dominant model '${DOM_MODEL}' does not exist or \
is not a file (option -D)." >&2
          exit 1
        fi
        ;; # -D <inFile>
      R)
        REC_MODEL="${OPTARG}"
        if [[ ! -f "${REC_MODEL}" ]]; then
          echo "[ERROR] Input recessive model '${REC_MODEL}' does not exist or \
is not a file (option -R)." >&2
          exit 1
        fi
        ;; # -R <inFile>

      m)
        MAP="${OPTARG}"
        if [[ ! -f "${MAP}" ]]; then
          echo "[ERROR] Input MAP file '${MAP}' does not exist or is not a file\
 (option -m)." >&2
          exit 1
        fi
        ;; # -m <inFile>
      d)
        DAT="${OPTARG}"
        if [[ ! -f "${DAT}" ]]; then
          echo "[ERROR] Input DAT file '${DAT}' does not exist or is not a file\
 (option -d)." >&2
          exit 1
        fi
        ;; # -d <inFile>
      p)
        PED="${OPTARG}"
        if [[ ! -f "${PED}" ]]; then
          echo "[ERROR] Input PED file '${PED}' does not exist or is not a file\
 (option -p)." >&2
          exit 1
        fi
        ;; # -p <inFile>
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
          echo "[ERROR] Output directory '${OUTDIR}' does not exist (option -o). \
Please create it." >&2
          exit 1
        fi
        ;; # -o <inDirectory>
      s)
        OUTPREFIX="${OPTARG}"
        ;; # -s <string>
      l)
        MINLODTH="${OPTARG}"
        if ! [[ "${MINLODTH}" =~ ^[0-9]+\.?[0-9]*$ ]] || \
          [[ "$(echo "${MINLODTH}<=0" | bc -l)" -eq 1 ]]; then
          echo '[ERROR] The lod-score threshold must be greater than 0 \
(option -l).' >&2
          exit 1
        fi
        ;; # -l <float>
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

  readonly DOM_MODEL REC_MODEL MAP DAT PED OUTDIR CHROMOSOME MINLODTH


  ### checking input directories and files
  if [[ "${DOM_MODEL}" = 'N.O.F.I.L.E' ]]; then
    echo '[ERROR] Input dominant model was not supplied (mandatory option -d)' >&2
    exit 1
  fi
  if [[ "${REC_MODEL}" = 'N.O.F.I.L.E' ]]; then
    echo '[ERROR] Input recessive model was not supplied (mandatory option -r)' >&2
    exit 1
  fi
  if [[ "${MAP}" = 'N.O.F.I.L.E' ]]; then
    echo '[ERROR] Input MAP file was not supplied (mandatory option -m)' >&2
    exit 1
  fi
  if [[ "${DAT}" = 'N.O.F.I.L.E' ]]; then
    echo '[ERROR] Input DAT file was not supplied (mandatory option -d)' >&2
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
    OUTPREFIX="results_sgl_chr"
  fi
  readonly OUTPREFIX


  ### print used parameters
  cat - <<EOF
    ${NAME}
    Parameters as interpreted:
      -D ${DOM_MODEL}
      -R ${REC_MODEL}
      -m ${MAP}
      -d ${DAT}
      -p ${PED}
      -c ${CHROMOSOME}
      -o ${OUTDIR}
      -s ${OUTPREFIX}
      -l ${MINLODTH}

EOF


  #main process
  ### output names
  OUTDOM="${OUTDIR}/${OUTPREFIX}${CHROMOSOME}_dominant.txt"
  OUTREC="${OUTDIR}/${OUTPREFIX}${CHROMOSOME}_recessive.txt"
  OUTDOMSIGNIF="${OUTDIR}/${OUTPREFIX}${CHROMOSOME}_dominant_LODsignif.txt"
  OUTRECSIGNIF="${OUTDIR}/${OUTPREFIX}${CHROMOSOME}_recessive_LODsignif.txt"


  ##run merlin if input files are not empty
  if [[ -s "${DAT}" ]] && [[ -s "${PED}" ]] && [[ -s "${MAP}" ]]; then
    echo "RUN SINGLEPOINT MERLIN ANALYSIS - CHROMOSOME ${CHROMOSOME}"
    merlin --quiet -d "${DAT}" -p "${PED}" -m "${MAP}" --model "${DOM_MODEL}" --singlepoint > "${OUTDOM}" \
      || {
        sed -n '/^FATAL ERROR/,$p' "${OUTDOM}" > "${OUTDOM}.err"
        echo "[ERROR] Merlin singlepoint analysis of the chromosome \
${CHROMOSOME} with a dominant model ends with error. See details into file \
${OUTDOM}.err" >&2
        END_WITH_ERRORS=1
    }
    merlin --quiet -d "${DAT}" -p "${PED}" -m "${MAP}" --model "${REC_MODEL}" --singlepoint > "${OUTREC}" \
      || {
        sed -n '/^FATAL ERROR/,$p' "${OUTREC}" > "${OUTREC}.err"
        echo "[ERROR] Merlin singlepoint analysis of the chromosome \
${CHROMOSOME} using a recessive model ends with error. See details into file \
${OUTREC}.err" >&2
        END_WITH_ERRORS=1
      }
  else
    echo "[WARNING] Ignoring merlin singlepoint analysis of the chromosome\
 ${CHROMOSOME}."
  fi;


  ##get significant LOD scores
  for model in dominant recessive; do
    file="${OUTDIR}/${OUTPREFIX}${CHROMOSOME}_${model}.txt"

    if [[ -f "${file}" ]] && [[ -s "${file}" ]]; then
      echo "[INFO] Get significant lod scores from merlin singlepoint analysis - \
chromosome ${CHROMOSOME}."
      prefix="${file%.txt}"
      line_nb="$(head -n -4 "${file}" | grep -n 'LOD ' | cut -f1 -d:)"
      if [[ ! -z "${line_nb}" ]]; then
        end_file="$(head -n -4 "$file" | wc -l | cut -f1 -d' ')"
        if [[ ! -z "${end_file}" ]]; then
          head -n -4 "${file}" \
            | tail -n "$(( ${end_file} - ${line_nb} ))" \
            > "${prefix}.woheader.txt"

          if [[ ! -f "${prefix}.woheader.txt" ]] \
            || [[ ! -s "${prefix}.woheader.txt" ]]; then
              echo "[ERROR] File ${prefix}.woheader.txt was not output or is empty." >&2
          fi

          if [[ -f "${prefix}.woheader.txt" ]] \
            && [[ -s "${prefix}.woheader.txt" ]]; then
              # print max lod score
              awk -v max=0 -v model="${model}" -v chr="${CHROMOSOME}" '{
                if($2>max){
                  max=$2
                }
              }END{
                print "chr"chr" - "model" model - singlepoint - max_observed_lod_score:\t"max
              }' "${prefix}.woheader.txt" \
                || echo "[ERROR] Program failed to display max observed lod score for \
chromosome ${CHROMOSOME}." >&2

              # get significant LOD scores
              awk -v th="${MINLODTH}" '{
                if($2>=th){print $0}
              }' "${prefix}.woheader.txt" \
                > "${prefix}_LODsignif.txt" \
                || echo '[ERROR] Program failed to retrieve significant lod scores.' >&2
          fi

          if [[ -f "${prefix}_LODsignif.txt" ]] \
            && [[ -s "${prefix}_LODsignif.txt" ]]; then
              echo "[INFO] Singlepoint analysis - ${model} model - chrom ${CHROMOSOME} \
- significant lod scores were output."
          fi

        fi
      else
        echo "[WARNING] Singlepoint analysis - ${model} model - chrom ${CHROMOSOME} \
- no lod score results were output."
        touch "${prefix}_LODsignif.txt"
      fi
    fi
  done


  #check output files
  if [[ ! -f "${OUTDOM}" ]] || [[ ! -s "${OUTDOM}" ]]; then
    echo "[ERROR] File '${OUTDOM}' was not output or is empty." >&2
  fi
  if [[ ! -f "${OUTREC}" ]] || [[ ! -s "${OUTREC}" ]]; then
    echo "[ERROR] File '${OUTREC}' was not output or is empty." >&2
  fi
  if [[ ! -f "${OUTDOMSIGNIF}" ]] || [[ ! -s "${OUTDOMSIGNIF}" ]]; then
    echo "[WARNING] File '${OUTDOMSIGNIF}' was not output or is empty."
  fi
  if [[ ! -f "${OUTRECSIGNIF}" ]] || [[ ! -s "${OUTRECSIGNIF}" ]]; then
    echo "[WARNING] File '${OUTRECSIGNIF}' was not output or is empty."
  fi

  if [[ "${END_WITH_ERRORS}" -eq 0 ]]; then
    echo 'DONE'
  else
    echo "[ERROR] chromosome ${CHROMOSOME} - some merlin singlepoint analyses \
failed" >&2
  fi

  return 0
}

main "$@"
