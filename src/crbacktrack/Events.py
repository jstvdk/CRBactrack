ev_TA_file = 'TA_gal_events.txt' # файл з ТА подіями
ev_PAO_file = 'PAO_events.txt'  # файл з Оже подіями


def extract_TA_events(ev_TA): # зчитування подій з файла, запис їх в 3 списки - енергія, довгота і широта, і формування результуючого списку списків
    TA_events = []
    file_w_events = open(ev_TA, 'r')  # open file with events
    events_list = file_w_events.readlines()  # read all lines into list of strings
    file_w_events.close()  # close events file
    E = []
    lon = []
    lat = []
    for i in events_list:
        sep_list = i.split()
        E.append(float(sep_list[0]))
        lon.append(float(sep_list[1]))
        lat.append(float(sep_list[2]))
    TA_events.append(E)
    TA_events.append(lon)
    TA_events.append(lat)
    return TA_events


def extract_PAO_events(ev_PAO): # зчитування подій з файла, запис їх в 3 списки - енергія, довгота і широта, і формування результуючого списку списків
    PAO_events = []
    file_w_events = open(ev_PAO, 'r')  # open file with events
    events_list = file_w_events.readlines()  # read all lines into list of strings
    file_w_events.close()  # close events file
    E = []
    lon = []
    lat = []
    for i in events_list:
        sep_list = i.split()
        E.append(float(sep_list[3]))
        lon.append(float(sep_list[6]))
        lat.append(float(sep_list[7]))
    PAO_events.append(E)
    PAO_events.append(lon)
    PAO_events.append(lat)
    return PAO_events


def all_events_no_calib(ta_events, pao_events):
    # TA * 1.0 and PAO * 1.0
    events_1020eV = []
    for i in range(len(ta_events[0])):
        if ta_events[0][i] > 100:
            events_1020eV.append([ta_events[0][i], ta_events[1][i], ta_events[2][i]])
    for i in range(len(pao_events[0])):
        if pao_events[0][i] > 100:
            events_1020eV.append([pao_events[0][i], pao_events[1][i], pao_events[2][i]])
    return events_1020eV


def all_events_ta1052_pao0948(ta_events, pao_events):
    # TA * 1.052 and PAO * 0.948
    events_1020eV = []
    for i in range(len(ta_events[0])):
        if ta_events[0][i] * 1.052 > 100:
            events_1020eV.append([ta_events[0][i] * 1.052, ta_events[1][i], ta_events[2][i]])
    for i in range(len(pao_events[0])):
        if pao_events[0][i] * 0.948 > 100:
            events_1020eV.append([pao_events[0][i] * 0.948, pao_events[1][i], pao_events[2][i]])
    return events_1020eV


def all_events_ta095_pao087(ta_events, pao_events):
    # TA / 0.95 and PAO / 0.87
    events_1020eV = []
    for i in range(len(ta_events[0])):
        if ta_events[0][i] / 0.95 > 100:
            events_1020eV.append([ta_events[0][i] / 0.95, ta_events[1][i], ta_events[2][i]])
    for i in range(len(pao_events[0])):
        if pao_events[0][i] / 0.87 > 100:
            events_1020eV.append([pao_events[0][i] / 0.87, pao_events[1][i], pao_events[2][i]])
    return events_1020eV


ta = extract_TA_events(ev_TA_file)
pao = extract_PAO_events(ev_PAO_file)

events_1020_no_calib = all_events_no_calib(ta, pao)
events_1020_ta1052_pao0948 = all_events_ta1052_pao0948(ta, pao)
events_1020_ta095_pao087 = all_events_ta095_pao087(ta, pao)

# for l in events_1020_no_calib:
#     print(" ".join(map(str, l)))
#
# print('\n')
#
# for l in events_1020_ta1052_pao0948:
#     print(" ".join(map(str, l)))
#
# print('\n')
#
# for l in events_1020_ta095_pao087:
#     print(" ".join(map(str, l)))