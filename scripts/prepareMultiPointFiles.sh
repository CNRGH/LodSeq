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
GENMAPS='N.O.D.I.R'               ## -g
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
    ${NAME} [options] -m <in_map> -p <in_ped> -c <string> -g <in_genetic_maps_dir> \
-o <out_dir> 
      -m <inFile>       input map file
      -p <inFile>       input ped pedigree file
      -c <string>       input chromosome
      -g <inDirectory>  input directory where are stored genetic map files \
without header
      -o <inDirectory>  directory where are stored output files
      -s <string>       prefix of output files
      -t <int>          number of threads used by plink steps (default : 1)
      -h                print help
  
  DESCRIPTION :
    ${NAME} formats .ped .dat and .map files for a multipoint merlin analysis from
    input .map .ped files.

  EXAMPLE :
    ${NAME} -m cohort_chr21.map -p cohort_chr21.ped -c 21 \
-g genetic_map_HapMapII_GRCh37/ -o outdir/ -s cohort_multi_chr -t 1
EOF

  return 0
}


################################################################################
# MAIN FUNCTION                                                                #
################################################################################
# main
# Formats .ped .dat and .map files for a multipoint merlin analysis from input 
#   .map .ped files.
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

  # if less than five arguments supplied, display usage
  if [[  "$#" -lt 5 ]]; then
    echo '[ERROR] Missing argument(s)' >&2
    display_usage
    exit 1
  fi



  ## catch option values
  while getopts :m:p:g:c:o:s:t: option
  do
    if [[ -z "${OPTARG}" ]]; then 
      echo "[ERROR] Empty argument for option -${option}" >&2 
      exit 1 
    fi

    case "${option}" in
      m)
        MAP="${OPTARG}"
        if [[ ! -f "${MAP}" ]]; then 
          echo "[ERROR] Input MAP file '${MAP}' does not exist (option -m)." >&2 
          exit 1  
        fi
        ;; # -m <inFile>
      p)
        PED="${OPTARG}"
        if [[ ! -f "${PED}" ]]; then 
          echo "[ERROR] Input PED '${PED}' does not exist (option -p)." >&2 
          exit 1  
        fi
        ;; # -p <inFile>
      g)
        GENMAPS="${OPTARG}"
        if [[ ! -d "${GENMAPS}" ]]; then 
          echo "[ERROR] Input directory of genetic maps '${GENMAPS}' does not \
exist (option -g)." >&2 
          exit 1  
        fi
        ;; # -g <inDirectory>
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

  readonly MAP PED GENMAPS OUTDIR CHROMOSOME THREADS


  ### checking input directories and files
  if [[ "${MAP}" = 'N.O.F.I.L.E' ]]; then 
    echo '[ERROR] Input MAP file was not supplied (mandatory option -m)' >&2 
    exit 1  
  fi
  if [[ "${PED}" = 'N.O.F.I.L.E' ]]; then 
    echo '[ERROR] Input PED file was not supplied (mandatory option -p)' >&2 
    exit 1  
  fi
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

  ### define default output file prefix
  if [[ "${OUTPREFIX}" = 'N.O.S.T.R.I.N.G' ]]; then 
    OUTPREFIX="$(basename "${MAP}" '.map')_multi"
  fi
  readonly OUTPREFIX


  ### print used parameters
  cat - <<EOF
    ${NAME} 
    Parameters as interpreted:
      -m ${MAP} 
      -p ${PED} 
      -c ${CHROMOSOME} 
      -g ${GENMAPS} 
      -o ${OUTDIR} 
      -s ${OUTPREFIX} 
      -t ${THREADS}

EOF


  #main process
  ##add genetic distances to map file
  echo 'PREPARE INPUT FILES OF MULTIPOINT MERLIN ANALYSIS'
 
  #because there is no genetic map for the Y chromosome
  if [[ "${CHROMOSOME}" != "Y" ]]; then 
    cp -p "${MAP}" "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map" || exit $? 
    cp -p "${PED}" "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped" || exit $?

    #first join to get the list of SNPs in common
    ##sort lexically not numerically, by position (bp), before join command
    awk '{ $4=sprintf("%015d", $4); print $0}' "${MAP}" \
      | sort -k 4 > "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.SORT.map" \
      || {
        echo '[ERROR] Sorting failed' >&2
        exit 1 
    }
    ## get numerically sorted list of SNPs - sort by position (bp)
    join -1 4 -2 1 "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.SORT.map" "${GENMAPS}/genetic_map_GRCh37_chr${CHROMOSOME}_wo_head.txt" \
      | sort -n -k 5 \
      | awk '{print $3}' \
      > "${OUTDIR}/list_snp_chr${CHROMOSOME}_with_gendist.txt" \
      || {
        echo "[ERROR] An error occurred when generating the list of SNPs in common \
'${OUTDIR}/list_snp_chr${CHROMOSOME}_with_gendist.txt', \
file does not exist or is empty." >&2 
        exit 1
    }
    rm "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.SORT.map"

    #check files exist
    if [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped" ]] \
      || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped" ]]; then 
        echo "[ERROR] File '${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped' does not \
exist or is empty." >&2 
        exit 1 
    fi
    if [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map" ]] \
      || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map" ]]; then 
        echo "[ERROR] File '${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map' does not \
exist or is empty." >&2 
        exit 1 
    fi
    if [[ ! -f "${OUTDIR}/list_snp_chr${CHROMOSOME}_with_gendist.txt" ]] \
      || [[ ! -s "${OUTDIR}/list_snp_chr${CHROMOSOME}_with_gendist.txt" ]]; then 
        echo "[ERROR] An error occurred when generating the list of SNPs in common \
'${OUTDIR}/list_snp_chr${CHROMOSOME}_with_gendist.txt', \
file does not exist or is empty." >&2 
        exit 1 
    fi

    #select only these SNPs in common in ped/map files
    plink_tmp_log="$(mktemp "${OUTDIR}/plink.XXX.log")"
    plink \
      --file "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK" \
      --extract "${OUTDIR}/list_snp_chr${CHROMOSOME}_with_gendist.txt" \
      --recode \
      --out "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK" \
      --threads "${THREADS}" \
      &> "${plink_tmp_log}" 
    rc=$?
    cat "${plink_tmp_log}" | tee >(sed -n '/^Error/,$p' >&2) | grep -v '^Error'
    if [[ -f "${plink_tmp_log}" ]]; then rm "${plink_tmp_log}"; fi
    if [[ $rc -ne 0 ]] \
      || [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.ped" ]] \
      || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.ped" ]] \
      || [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.map" ]] \
      || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.map" ]]; then
        echo "[ERROR] plink failed while generating the files \
'${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.ped' \
'${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.map'. \
See log files for more details." >&2
        exit $rc
    fi 
      
      
    #second join to add the genetic distances to the .map files
    ##sort lexically not numerically, before join command
    awk '{ 
      $4=sprintf("%015d", $4); 
      print $0
    }' "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.map" \
      | sort -k 4 \
      > "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.TMPSORT.map" \
      || {
        echo '[ERROR] Sorting failed' >&2
        exit 1
    }

    ## get numerically sorted .map and .dat files as inputs of merlin
    join -1 4 -2 1 "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.TMPSORT.map" "${GENMAPS}/genetic_map_GRCh37_chr${CHROMOSOME}_wo_head.txt" \
      | sort -n -k 5 \
      | awk 'BEGIN{OFS="\t"};{print $2,$3,$5}' \
      > "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.JOIN.map" \
      || {
        echo "[ERROR] An error occurred when generating the map file \
'${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.JOIN.map'" >&2 
        exit 1
    }

    # remove tmp files 
    rm "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.TMPSORT.map"
    rm "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.map"
    rm "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.ped"
    rm "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.INPLINK.map"

    #make a .dat file for merlin with bi-allelic markers with a genetic distance
    echo 'A affection' > "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat" \
      && awk '{print "M",$2}' "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.JOIN.map" \
      >> "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat" \
      || {
        echo "[ERROR] An error occurred when generating the dat file \
${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat'" >&2 
        exit 1
    }

    mv "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.JOIN.map" "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map" || exit $?
    mv "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.ped" "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.ped" || exit $?

    #rename plink log file
    if [[ -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.log" ]]; then 
      mv "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.OUTPLINK.log" "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.log" 
    fi

  elif [[ "${CHROMOSOME}" == "Y" ]]; then
    echo "[INFO] Ignoring merlin multipoint analysis of the chromosome ${CHROMOSOME}."
    touch "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat"
    touch "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map"
    touch "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.ped"
  fi

  #remove temp files
  if [[ -f "${OUTDIR}/list_snp_chr${CHROMOSOME}_with_gendist.txt" ]]; then 
    rm "${OUTDIR}/list_snp_chr${CHROMOSOME}_with_gendist.txt" 
  fi

  #check output files
  if [[ "${CHROMOSOME}" != "Y" ]]; then
    if [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map" ]] \
      || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map" ]]; then
        echo "[ERROR] File '${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.map' was not output \
or is empty." >&2
        exit 1
    fi
    if [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.ped" ]] \
      || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.ped" ]]; then
        echo "[ERROR] File '${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.ped' was not output \
or is empty." >&2
        exit 1
    fi
    if [[ ! -f "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat" ]] \
      || [[ ! -s "${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat" ]]; then
        echo "[ERROR] File '${OUTDIR}/${OUTPREFIX}${CHROMOSOME}.dat' was not output \
or is empty." >&2
        exit 1
    fi
  fi
 
  return 0  
}

main "$@"
