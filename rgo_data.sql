select distinct substr(CSG,1,4) CSG,
substr(date,9,2) day, 
truncate(round((CWSA),0),0) CWSA,
truncate(round ((LNS),0),0) LNS,
truncate(round ((CLD),0),0) CLD, 
truncate(round((CUA),0),0) CUA INTO OUTFILE 'c:/sun/rgo_data_5809.csv' FIELDS TERMINATED BY ','
from rgo_data 
where substr(date,1,4) >= '1958'
and substr(date,1,4) <= '1958'
and substr(date,7,1) >= '9'
and substr(date,7,1) <= '9'
;

