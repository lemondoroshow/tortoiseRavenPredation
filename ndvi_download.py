import os
import pandas as pd
import selenium
import time
from datetime import datetime as dt
from pathlib import Path
from selenium import webdriver
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait


def wait_for_download(dl_path):
    wait = True
    while wait:
        time.sleep(1)
        wait = False
        for file in os.listdir(dl_path):
            if file.endswith('.crdownload'):
                wait = True

# Find which NDVIs are already downloaded and processed
downloaded = os.listdir(Path.cwd() / "data/ndvi/processed/")
downloaded = sorted([int(d[:-4].replace("-", "")) for d in downloaded])

# Get list of dates including and after 2001
all_dates = pd.read_csv("data/bbs_obs_mojave.csv")
dates = [int(dt.strftime(dt.strptime(d, "%Y-%m-%d"), format = "%Y%m%d"))
             for d in all_dates["Date"]]
dates = sorted(list(set([d for d in dates if 
                        (d >= 20010101 and d not in downloaded)])))

# Iterate through years
print("Enter start year: ")
start = int(input())
print("Enter end year: ")
end = int(input()) + 1
for year in range(start, end):
    doi = ["_" + str(d) + "_" for d in dates if 
           (d >= year * 10**4 and d < (year + 1) * 10**4)]

    # Set webdriver settings
    download_folder = Path.cwd() / ("data/ndvi/raw/" + str(year) + "/")
    chr_options = Options()
    prefs = {"download.default_directory" : str(download_folder)}
    chr_options.add_experimental_option("prefs", prefs)
    chr_options.add_argument('--headless=new')
    driver = webdriver.Chrome(options = chr_options)

    # Iterate through all dates of interest
    for date in doi:

        # Debugging
        print(date)

        # Reset driver
        driver.get("https://noaa-cdr-ndvi-pds.s3.amazonaws.com/index.html#data/" + str(year) + "/")
        while len(driver.find_elements(By.ID, "dt-search-0")) == 0:
            WebDriverWait(1, timeout = 5)

        # Filter for date
        search_box = driver.find_element(By.ID, "dt-search-0")
        search_box.send_keys(date)

        # Wait until download is there
        while len(driver.find_elements(By.XPATH, "/html/body/div/div/div/div/div[2]/div/div[2]/div/table/tbody/tr[1]/td[1]/a")) == 0:
            WebDriverWait(1, timeout = 5)

        # Get download link
        link = driver.find_element(By.XPATH, "/html/body/div/div/div/div/div[2]/div/div[2]/div/table/tbody/tr[1]/td[1]/a")
        link.click()

        # Wait for download to finish
        wait_for_download(download_folder)

        # Clear search box
        for i in range(0, 11):
            search_box.send_keys(Keys.BACKSPACE)

    # Quit driver
    driver.quit()
