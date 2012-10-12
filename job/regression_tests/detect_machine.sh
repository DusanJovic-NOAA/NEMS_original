case `hostid` in

  0xac7c006)          MACHINE_ID=ccs ;;     ### cirrus1
  0xac7c012)          MACHINE_ID=ccs ;;     ### cirrus2

  0xac7a006)          MACHINE_ID=ccs ;;     ### stratus1
  0xac7a012)          MACHINE_ID=ccs ;;     ### stratus2

  0xffffffffbcc090b2) MACHINE_ID=gaea ;;    ### gaea1
  0xffffffffbcc08fb2) MACHINE_ID=gaea ;;    ### gaea2
  0xffffffffbcc08eb2) MACHINE_ID=gaea ;;    ### gaea3
  0xffffffffbcc08db2) MACHINE_ID=gaea ;;    ### gaea4

  ae0a2500)           MACHINE_ID=zeus ;;    ### zeus1
  ae0a2600)           MACHINE_ID=zeus ;;    ### zeus2
  ae0a2700)           MACHINE_ID=zeus ;;    ### zeus3
  ae0a2800)           MACHINE_ID=zeus ;;    ### zeus4
  ae0a2900)           MACHINE_ID=zeus ;;    ### zeus5
  ae0a2a00)           MACHINE_ID=zeus ;;    ### zeus6
  ae0a7d13)           MACHINE_ID=zeus ;;    ### zeus7
  ae0a7e13)           MACHINE_ID=zeus ;;    ### zeus8

  b20a057e)           MACHINE_ID=jet ;;     ### jet1
  b20a067e)           MACHINE_ID=jet ;;     ### jet2
  b20a077e)           MACHINE_ID=jet ;;     ### jet3
  b20a087e)           MACHINE_ID=jet ;;     ### jet4

  010afa01)           MACHINE_ID=eddy ;;    ### eddy

esac


# --- for Zeus, find available account ID
if [ ${MACHINE_ID} = zeus ]; then
  clear
  ACCNR=null
  for i in rm gm ada cmp gefs sepp omd
  do
    printf %s "Looking for accout ID: " $i " ..."
    nr=`account_params 2>&1 | grep -v Initial | grep " "$i" " | wc -l`

    if [ $nr -eq 1 ]; then
      ACCNR=$i ; echo OK
    else
      echo
    fi
  done
  if [ $ACCNR = null ]; then echo "Check your account ID"; exit ; fi
  clear
fi
