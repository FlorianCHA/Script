if record.CHROM in dicoRNA.keys() and record.POS in dicoRNA[record.CHROM].keys():
    listeId = dicoRNA[record.CHROM][record.POS].keys()

    if record.POS in listePosition :
        for id in listeId :
            if len(record.REF) > len(ac):
                ac_new = recupAC(id, record, MQmin, DPmin, QDmin)
                dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac_new)
            else:
                dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac)
        for i in range(1, len(record.REF)):
                listePosition.append(record.POS + i)
    else :
        listePosition = [record.POS]
        for id in listeId :
            brin = dicoRNA[record.CHROM][record.POS][id][1]
            if record.INFO == {}:
                ac = 'N'
                if id not in dicoSeq.keys():
                    dicoSeq[id] = ac
                    dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
                else:
                    dicoSeq[id] += ac
            else:
                brin = dicoRNA[record.CHROM][record.POS][id][1]
                if id not in dicoSeq.keys():
                    ac = recupAC(id, record, MQmin, DPmin, QDmin)
                    dicoSeq[id] = ac
                    dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
                else:
                    ac = recupAC(id, record, MQmin, DPmin, QDmin)
                    dicoSeq[id] += ac





    for id in listeId:
        # listePosition = [0] # ATTENTION PROBLEME DE listePOSITION réinitialisé comme avant
        if record.POS in listePosition and id in dicoSeq.keys():
            if len(record.REF) > len(ac):
                ac_new = recupAC(id, record, MQmin, DPmin, QDmin)
                dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac_new)

            else:
                dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac)
            for i in range(1, len(record.REF)):
                listePosition.append(record.POS + i)

        # elif record.POS > int(listePosition[-1]) + 1 and id in dicoSeq.keys() :
        #     dicoSeq[id] += 'N' * ((int(listePosition[-1]) + 1) - int(record.POS))
        #
        # elif id not in dicoSeq.keys() and int(dicoRNA[record.CHROM][record.POS][id][3]) < record.POS :
        #     dicoSeq[id] = 'N' * (record.POS - int(record.POS - dicoRNA[record.CHROM][record.POS][id][3]))

        else:
            listePosition = [record.POS]
            brin = dicoRNA[record.CHROM][record.POS][id][1]
            if record.INFO == {}:
                ac = 'N'
                if id not in dicoSeq.keys():
                    dicoSeq[id] = ac
                    dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
                else:
                    dicoSeq[id] += ac
            else:
                brin = dicoRNA[record.CHROM][record.POS][id][1]
                if id not in dicoSeq.keys():
                    ac = recupAC(id, record, MQmin, DPmin, QDmin)
                    dicoSeq[id] = ac
                    dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
                else:
                    ac = recupAC(id, record, MQmin, DPmin, QDmin)
                    dicoSeq[id] += ac










  if record.CHROM in dicoRNA.keys() and record.POS in dicoRNA[record.CHROM].keys():
                listeId = dicoRNA[record.CHROM][record.POS].keys()
                for id in listeId:
                    # listePosition = [0] # ATTENTION PROBLEME DE listePOSITION réinitialisé comme avant
                    if record.POS in listePosition and id in dicoSeq.keys():
                        if len(record.REF) > len(ac):
                            ac_new = recupAC(id, record, MQmin, DPmin, QDmin)
                            dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac_new)

                        else:
                            dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac)
                        for i in range(1, len(record.REF)):
                            listePosition.append(record.POS + i)

                    # elif record.POS > int(listePosition[-1]) + 1 and id in dicoSeq.keys() :
                    #     dicoSeq[id] += 'N' * ((int(listePosition[-1]) + 1) - int(record.POS))
                    #
                    # elif id not in dicoSeq.keys() and int(dicoRNA[record.CHROM][record.POS][id][3]) < record.POS :
                    #     dicoSeq[id] = 'N' * (record.POS - int(record.POS - dicoRNA[record.CHROM][record.POS][id][3]))

                    else:
                        listePosition = [record.POS]
                        brin = dicoRNA[record.CHROM][record.POS][id][1]
                        if record.INFO == {}:
                            ac = 'N'
                            if id not in dicoSeq.keys():
                                dicoSeq[id] = ac
                                dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
                            else:
                                dicoSeq[id] += ac
                        else:
                            brin = dicoRNA[record.CHROM][record.POS][id][1]
                            if id not in dicoSeq.keys():
                                ac = recupAC(id, record, MQmin, DPmin, QDmin)
                                dicoSeq[id] = ac
                                dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
                            else:
                                ac = recupAC(id, record, MQmin, DPmin, QDmin)
                                dicoSeq[id] += ac







































if record.CHROM in dicoRNA.keys() and record.POS in dicoRNA[record.CHROM].keys():
    listeId = dicoRNA[record.CHROM][record.POS].keys()
    for id in listeId :
        if record.POS in listePosition:
            if len(record.REF) > len(ac):
                ac_new = recupAC(id, record, MQmin, DPmin, QDmin)
                dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac_new)

            else:
                dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac)
            for i in range(1, len(record.REF)):
                listePosition.append(record.POS + i)
        elif record.POS > int(listePosition[-1]) + 1 and id in dicoSeq.keys() :
             dicoSeq[id] = 'N' * ((int(listePosition[-1]) + 1) - int(record.POS)
        elif id not in dicoSeq.keys() and int(dicoRNA[record.CHROM][record.POS][id][3]) < record.POS :
            dicoSeq[id] = 'N'*(record.POS - int(record.POS - dicoRNA[record.CHROM][record.POS][id][3]))


        else:
             listePosition =[record.POS]
             brin = dicoRNA[record.CHROM][record.POS][id][1]
            if record.INFO == {}:
                ac = 'N'

            if id not in dicoSeq.keys():
                dicoSeq[id] = ac
                dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
            else:
                dicoSeq[id] += ac

            else:
            brin = dicoRNA[record.CHROM][record.POS][id][1]
            if id not in dicoSeq.keys():
                ac = recupAC(id, record, MQmin, DPmin, QDmin)
                dicoSeq[id] = ac
                dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
            else:
                ac = recupAC(id, record, MQmin, DPmin, QDmin)
                dicoSeq[id] += ac



             ----------------------
             if record.CHROM in dicoRNA.keys() and record.POS in dicoRNA[record.CHROM].keys():
                 if
             record.POS in listePosition:
             if len(record.REF) > len(ac):
                 ac_new = recupAC(id, record, MQmin, DPmin, QDmin)
             for id in listeId:
                 dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac_new)

             else:
                 for
             id in listeId:
             dicoSeq[id] = dicoSeq[id][:-len(ac)] + 'N' * len(ac)
             for i in range(1, len(record.REF)):
                 listePosition.append(record.POS + i)
             elif record.POS > listePosition[-1] + 1:
             for id in listeID:
                 if
             id in dicoSeq.keys():
             dicoSeq[id] = 'N' * ((int(listePosition[-1]) + 1) - int(record.POS)




                                  else:
                                  listePosition =[record.POS]
                                  listeId = dicoRNA[record.CHROM][record.POS].keys()

             if record.INFO == {}:
                 ac = 'N'
             for id in listeId:
                 brin = dicoRNA[record.CHROM][record.POS][id][1]
             if id not in dicoSeq.keys():
                 dicoSeq[id] = ac
             dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
             else:
             dicoSeq[id] += ac

             else:
             for id in listeId:
                 brin = dicoRNA[record.CHROM][record.POS][id][1]
             if id not in dicoSeq.keys():
                 ac = recupAC(id, record, MQmin, DPmin, QDmin)
             dicoSeq[id] = ac
             dicoInfo[id] = dicoRNA[record.CHROM][record.POS][id]
             else:
             ac = recupAC(id, record, MQmin, DPmin, QDmin)
             dicoSeq[id] += ac