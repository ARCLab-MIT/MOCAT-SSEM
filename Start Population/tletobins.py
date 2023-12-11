# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 09:43:31 2022

@author: miles
"""
import pandas as pd
import numpy as np
import requests
from datetime import datetime
import os
import matplotlib.pyplot as plt
from dotenv import load_dotenv, find_dotenv

dir_path = os.path.dirname(os.path.realpath(__file__))
dotenv_path = os.path.join(dir_path, ".env")
if load_dotenv(dotenv_path=dotenv_path):
    account = os.environ.get("email")
    password = os.environ.get("password")
else:
    print("Could not find a .env file."
          " Please create a file called .env with your space-track.org email"
          " and password.")
    print("\nemail=your@email.com"
          "\npassword=yourPassword")

pd.set_option('display.max_rows', 50)
pd.set_option('display.max_columns', 50)
pd.set_option('display.width', 1000)
s = requests.Session()

def get_cookie(account, password):
  data = {'identity':account,'password':password}
  login_url = 'https://www.space-track.org/ajaxauth/login'
  response = s.post(login_url, data = data)
  login_cookie = response.cookies.get_dict()
  print("Get Cookie Response " + str(response))
  return login_cookie

def getTLES(login_cookie, limit = 100, apoapsis_limit = 2000):
    """
    """
    limit = str(limit)
    baseLink = 'https://www.space-track.org/basicspacedata/query/class/gp/'
    #decay = 'decay_date' + "/" + "null-val"
    epoch = 'EPOCH/' + "%3Enow-30" "/"
    alt = '/APOAPSIS/<{apoapsis_limit}/'.format(apoapsis_limit=apoapsis_limit)
    other = 'orderby/norad_cat_id/limit/' + limit + '/format/csv/emptyresult/show'
       
    link = baseLink  + epoch + alt + other
    print(link)
    r = s.get(link, cookies = login_cookie)
    savePath = "temp_GP_data_" + datetime.now().isoformat().replace(":","_") + ".csv"
    with open(savePath, "wb") as f:
        f.write(r.content)
    df = pd.read_csv(savePath)
        
    return df

def getTLEBins(bin_width, start_alt, end_alt, savedir = None, graph = False):
    # bin_width = altitude shell width in km
    # start_alt = starting altitude above the Earth in km
    # end_alt = ending altitude above the Earth in km
    # savedir = string with the folder where to save results files.
    # graph = bool for printing graphs or not.
    
    if savedir is not None:
        os.chdir(savedir)
    if savedir is None:
        savedir = os.getcwd()
        
    r_Earth = 6378.1366 #km
    login_cookie = get_cookie(account, password)
    data = getTLES(login_cookie, limit = 'null-val')
    #filterData =  data[[ 'OBJECT_NAME', 'MEAN_ELEMENT_THEORY', 'EPOCH', 'MEAN_MOTION', 'ECCENTRICITY', 'INCLINATION', 'RA_OF_ASC_NODE', 'MEAN_ANOMALY', 'NORAD_CAT_ID', 'SEMIMAJOR_AXIS','OBJECT_TYPE']]

    dists = []
    for objectType in data.OBJECT_TYPE.unique():
        bins = pd.interval_range(start=r_Earth + start_alt, end=r_Earth + end_alt, freq = bin_width)
        tempData = data[data.OBJECT_TYPE == objectType].copy(deep = True) 
        if objectType == 'PAYLOAD':
            constellations = [r"STARLINK", r"ORBCOMM", r"IRIDIUM", r"ONEWEB"]  
            masks = []
            for constellation in constellations:
                print(constellation)
                mask = tempData.OBJECT_NAME.str.contains(constellation)
                masks.append(mask)
            tempDataSlot = tempData[np.logical_or.reduce(masks)]
            tempDataUnSlot = tempData[~np.logical_or.reduce(masks)]
            
            tempDataSlot.loc[:,'range'] = pd.cut(tempDataSlot.SEMIMAJOR_AXIS,bins = bins).copy()
            tempDataSlot.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count').to_dict(orient='records')
            counts = tempDataSlot.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count')
            counts["range"] = counts["range"].astype("object")
            counts.to_csv("Counts_" + objectType + "slot" + "_bins_" + str(bin_width) + ".csv", index = False)
            dists.append(counts)

            tempDataUnSlot.loc[:,'range'] = pd.cut(tempDataUnSlot['SEMIMAJOR_AXIS'],bins = bins).copy(deep=True)
            tempDataUnSlot.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count').to_dict(orient='records')
            counts = tempDataUnSlot.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count')
            counts["range"] = counts["range"].astype("object")
            counts.to_csv("Counts_" + objectType + "unslot" + "_bins_" + str(bin_width) + ".csv", index = False)
            dists.append(counts)

        elif objectType == 'DEBRIS':  
            tempDataDerelict = tempData[(tempData.RCS_SIZE=="LARGE")]
            #tempDataDerelict = tempData[(tempData.RCS_SIZE=="LARGE") | (tempData.RCS_SIZE=="MEDIUM")]
            tempDataDerelict.loc[:,'range'] = pd.cut(tempDataDerelict.SEMIMAJOR_AXIS,bins = bins).copy(deep=True)
            tempDataDerelict.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count').to_dict(orient='records')
            counts = tempDataDerelict.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count')
            counts["range"] = counts["range"].astype("object")
            counts.to_csv("Counts_" + "DERELICT" + "_bins_" + str(bin_width) + ".csv", index = False)
            dists.append(counts)    

            tempDataDebris= tempData[(tempData.RCS_SIZE!="LARGE")]            
            #tempDataDebris= tempData[(tempData.RCS_SIZE!="LARGE") & tempData.RCS_SIZE!="MEDIUM"]
            tempDataDebris.loc[:,'range'] = pd.cut(tempDataDebris.SEMIMAJOR_AXIS,bins = bins).copy(deep=True)
            tempDataDebris.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count').to_dict(orient='records')
            counts = tempDataDebris.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count')
            counts["range"] = counts["range"].astype("object")
            counts.to_csv("Counts_" + objectType + "_bins_" + str(bin_width) + ".csv", index = False)
            dists.append(counts)    
        else:
            tempData.loc[:,'range'] = pd.cut(tempData.SEMIMAJOR_AXIS,bins = bins).copy(deep=True)
            tempData.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count').to_dict(orient='records')
            counts = tempData.groupby('range')['SEMIMAJOR_AXIS'].count().reset_index(name='Count')
            counts["range"] = counts["range"].astype("object")
            counts.to_csv("Counts_" + objectType + "_bins_" + str(bin_width) + ".csv", index = False)
            dists.append(counts)
            
    if graph == True:
        fig, ax = plt.subplots()
        labels = ['S','$S_U$', 'ROCKET BODY', 'DERELICT', 'DEBRIS'] #, 'UNKNOWN']
        for count, label in zip(dists, labels):
            #dists = dists.astype("interval")
            count["Avg"] = [(count.left + count.right)/2 for count in counts["range"].astype("interval")]
            ax.plot(count["Avg"] - r_Earth, count.Count, label = label )
        ax.set_xlabel("Altitude Bin Center [km]")
        ax.set_ylabel("Bin Population Count")
        #ax.set_ylim(0,100*bin_width)
        ax.legend()
        today = datetime.today().isoformat().replace(":","_")
        figname = "All_width_{bin_width}km_start_{start_alt}km_end_{end_alt}km_{today}.png".format(bin_width=bin_width, start_alt=start_alt, end_alt = end_alt, today=today)
        fig.savefig(fname=figname)
        plt.show()
        
        labels = ['S','$S_U$', 'ROCKET BODY', 'DERELICT', 'DEBRIS']#, 'UNKNOWN']
        for count, label in zip(dists, labels):
            fig, ax = plt.subplots()
            #dists = dists.astype("interval")
            count["Avg"] = [(count.left + count.right)/2 for count in counts["range"].astype("interval")]
            ax.plot(count["Avg"] - r_Earth, count.Count, label = label )
            ax.legend()
            ax.set_xlabel("Altitude Bin Center [km]")
            ax.set_ylabel("Bin Population Count")
            label = label.replace("$","") #Get rid of $ in file name
            figname = "{label}_width_{bin_width}km_start_{start_alt}km_end_{end_alt}km_{today}.png".format(label=label, bin_width=bin_width, start_alt=start_alt, end_alt = end_alt, today = today)
            fig.savefig(fname= figname)
            plt.show()
        return savedir
#%%

if __name__ == "__main__":
    getTLEBins(bin_width=50, start_alt=200, end_alt=2000, graph = True)