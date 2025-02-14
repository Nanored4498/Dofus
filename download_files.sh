#!/bin/env bash

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR

wget https://raw.githubusercontent.com/Trameurs/DofusFashionistaVanced/refs/heads/main/itemscraper/all_sets_fr.json
wget https://github.com/Trameurs/DofusFashionistaVanced/raw/refs/heads/main/itemscraper/all_equipment_fr.json 
wget https://raw.githubusercontent.com/Trameurs/DofusFashionistaVanced/refs/heads/main/itemscraper/all_mounts_fr.json