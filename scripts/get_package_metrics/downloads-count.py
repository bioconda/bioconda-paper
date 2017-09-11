#!/usr/bin/env python

import os, time, csv
import json
import dateutil.parser
from datetime import datetime, timezone
import urllib.request
import concurrent.futures
from collections import OrderedDict


def get_package_metrics(package):
    package_info = {}

    package_date = datetime.now(timezone.utc) #datetime.utcfromtimestamp(0)
    package_newest_build = dateutil.parser.parse("1970-01-01 00:00:00.000000+00:00")
    has_linux_build = False
    has_osx_build   = False

    # the json request includes a url that is sometimes wrong...
    preq = urllib.request.Request(package["url"].replace("anaconda.org/packages","anaconda.org/package"))
    #time.sleep(0.25) # to avoid being rate limited
    with urllib.request.urlopen(preq) as presponse:
        presult = json.loads(presponse.read().decode('utf-8'))
        total_downloads = 0
        for build in presult["files"]:
            if include_all_versions or (build["version"] == package["versions"][-1]):
                total_downloads += int(build["ndownloads"])

            build_date = dateutil.parser.parse(build["upload_time"])

            # set the package date if it is older than the oldest seen so far
            if build_date < package_date:
                package_date = build_date

            if build_date > package_newest_build:
                package_newest_build = build_date

            if build["attrs"]["operatingsystem"] == "linux" and not has_linux_build:
                has_linux_build = True

            if build["attrs"]["operatingsystem"] == "osx" and not has_osx_build:
                has_osx_build = True
            

        timedelta_available = (datetime.now(timezone.utc) - package_date)
        days_available      = timedelta_available.days + ((timedelta_available.seconds+60)/(3600*24)) # add a minute to avoid div by zero
        
        timedelta_to_newest = (datetime.now(timezone.utc) - package_newest_build) 
        days_since_newest   = timedelta_to_newest.days + ((timedelta_to_newest.seconds+60)/(3600*24))

        print("total downloads for",package["name"]+":",total_downloads)
        print("total days_available for",package["name"]+":",days_available)

        downloads_per_day = total_downloads / days_available
        print("total downloads_per_day for",package["name"]+":",downloads_per_day)
        print("")

        package_info["name"]              = package["name"]
        package_info["total_downloads"]   = total_downloads
        package_info["days_available"]    = days_available
        package_info["downloads_per_day"] = downloads_per_day
        package_info["days_since_newest"] = days_since_newest
        package_info["has_linux_build"]   = 1 if has_linux_build else 0
        package_info["has_osx_build"]     = 1 if has_osx_build else 0

        return package_info

req = urllib.request.Request('https://api.anaconda.org/packages/bioconda')

include_all_versions = True

total_downloads_all = 0
repo_info = {}

with urllib.request.urlopen(req) as response:
    result = json.loads(response.read().decode('utf-8'))

    with concurrent.futures.ProcessPoolExecutor() as executor:
        for execution in executor.map(get_package_metrics, result, chunksize=50):
                package_info = execution
                total_downloads_all += package_info["total_downloads"]
                repo_info[package_info["name"]] = package_info

with open("download_counts.csv", "wt") as outf:
    fieldnames = "name,total_downloads,days_available,downloads_per_day,days_since_newest,has_linux_build,has_osx_build".split(",")

    writer = csv.DictWriter(outf, fieldnames=fieldnames)

    writer.writeheader()    

    for package, package_info in OrderedDict(sorted(repo_info.items(), key=lambda t: t[1]["downloads_per_day"], reverse=True)).items():
        writer.writerow(package_info)

print("Total number of packages:", len(repo_info.keys()))
print("Total number of downloads:", total_downloads_all)

