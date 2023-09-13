#!/bin/bash

MD_RUN=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

cat sugars_list.txt | while read sugars

do

if [[ "$sugars" == *\/* ]]; then


cat bases_list.txt | while read ru

do

echo "copying $ru mol file to $sugars/$ru location"

cp $MD_RUN/vacuum_files/$ru.mol $MD_RUN/vacuum/$sugars/$ru/

cp $MD_RUN/water_files/$ru.mol $MD_RUN/water/$sugars/$ru/


cd ../../
done

else

cat bases_list.txt | while read ru

do

echo "copying $ru mol file to $sugars/$ru location"

cp $MD_RUN/vacuum_files/$ru.mol $MD_RUN/vacuum/$sugars/$ru/

cp $MD_RUN/water_files/$ru.mol $MD_RUN/water/$sugars/$ru/

cd ../
done
fi
done









