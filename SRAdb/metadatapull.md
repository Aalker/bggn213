---
title: "Getting Latitude and Longitude with Reutils"
author: "Amanda Alker"
date: "6/19/2018"
output: 
  html_document: 
    keep_md: yes
---




##Start with getting the right IDs and Search Queries



```r
#load the package
library(rentrez)
#This allows you to check out the databases and fields that you want to search
library(reutils)
#This is where the functions for getting db summaries and crosslinking the IDs are found
library(XML)
library(xml2)
#These will help you parse the output
```


Check out the Searchable space...



```r
#Gives the list of searchable databases
entrez_dbs()
```

```
##  [1] "pubmed"          "protein"         "nuccore"        
##  [4] "ipg"             "nucleotide"      "nucgss"         
##  [7] "nucest"          "structure"       "sparcle"        
## [10] "genome"          "annotinfo"       "assembly"       
## [13] "bioproject"      "biosample"       "blastdbinfo"    
## [16] "books"           "cdd"             "clinvar"        
## [19] "clone"           "gap"             "gapplus"        
## [22] "grasp"           "dbvar"           "gene"           
## [25] "gds"             "geoprofiles"     "homologene"     
## [28] "medgen"          "mesh"            "ncbisearch"     
## [31] "nlmcatalog"      "omim"            "orgtrack"       
## [34] "pmc"             "popset"          "probe"          
## [37] "proteinclusters" "pcassay"         "biosystems"     
## [40] "pccompound"      "pcsubstance"     "pubmedhealth"   
## [43] "seqannot"        "snp"             "sra"            
## [46] "taxonomy"        "biocollections"  "unigene"        
## [49] "gencoll"         "gtr"
```



These are the searchable fields within the database of interest. If anything, perform the search with one example on the web interface to see where you need to be looking.




```r
#Gives the list of searchable terms within database 'sra'

entrez_db_searchable("sra")
```

```
## Searchable fields for database 'sra'
##   ALL 	 All terms from all searchable fields 
##   UID 	 Unique number assigned to publication 
##   FILT 	 Limits the records 
##   ACCN 	 Accession number of sequence 
##   TITL 	 Words in definition line 
##   PROP 	 Classification by source qualifiers and molecule type 
##   WORD 	 Free text associated with record 
##   ORGN 	 Scientific and common names of organism, and all higher levels of taxonomy 
##   AUTH 	 Author(s) of publication 
##   PDAT 	 Date sequence added to GenBank 
##   MDAT 	 Date of last update 
##   GPRJ 	 BioProject 
##   BSPL 	 BioSample 
##   PLAT 	 Platform 
##   STRA 	 Strategy 
##   SRC 	 Source 
##   SEL 	 Selection 
##   LAY 	 Layout 
##   RLEN 	 Percent of aligned reads 
##   ACS 	 Access is public or controlled 
##   ALN 	 Percent of aligned reads 
##   MBS 	 Size in megabases
```



Based on the Ids we have, it looks like we will have to find 'Biosample' Ids from the SRA run Ids. To do that, we will need to look up the uids for each database and convert them. 


#Convert the Ids


First, let's look for the uids that we want. This is just a practice with one.




```r
#make the 'uid' object to find the specific run ids in the sra database

uid <- esearch("SRR3719567", db = "sra")
uid
```

```
## Object of class 'esearch' 
## List of UIDs from the 'sra' database.
## [1] "2684569"
```

```r
#careful! These search queries are hard-coded, so we need to figure out a better way to have inputs for the function
```



So it looks like we got something out as a uid. Go and check it online to make sure you are getting the correct uid out..... it works!




```r
#find 'sra uid' from the SRA run id
sra_uid <- elink(uid, dbFrom = "sra", dbTo = "biosample" )
sra_uid
```

```
## Object of class 'elink' 
## ELink query from database 'sra' to destination database 'biosample'.
## Query UIDs:
## [1] "2684569"
## Summary of LinkSet:
##        DbTo      LinkName LinkCount
## 1 biosample sra_biosample         1
```



So now lets convert the sra_uid to the biosample_uid. Then let's double check that the right biosample_uid is coming up too. 





```r
#Reveal the biosample_uid that will be used to access the attribute metadata
biosample_uid <- linkset(sra_uid)
biosample_uid
```

```
## $sra_biosample
## List of linked UIDs from database 'sra' to 'biosample'.
## [1] "5224525"
```



Great! Now let's try to find the summary page of the biosample uid


#Fetch the summary data 




```r
biosample_data <- esummary(biosample_uid, db = "biosample")
```



ITS IN THERE! Now we just have to parse it out!



##Parsing the data out of XML




```r
parsed_biosample_full <- content(biosample_data, "parsed")
View(parsed_biosample_full)
```



Okay, we're on to something, but it looks like the SampleData column did not parse out all of the way. Let's keep on working on it.




```r
parsed_biosample_sampledata <- xmlParseString(parsed_biosample_full$SampleData)
parsed_biosample_sampledata
```

```
## <BioSample access="public" publication_date="2016-06-09T17:54:29.673" last_update="2017-07-12T15:57:48.613" submission_date="2016-06-09T17:54:30.157" id="5224525" accession="SAMN05224525">
##   <Ids>
##     <Id db="BioSample" is_primary="1">SAMN05224525</Id>
##     <Id db="SRA">SRS1525469</Id>
##     <Id db="DOE Joint Genome Institute">Gp0091327</Id>
##   </Ids>
##   <Description>
##     <Title>Marine microbial communities from expanding oxygen minimum zones in the Saanich Inlet - SI073_LV_120m_DNA</Title>
##     <Organism taxonomy_id="408172" taxonomy_name="marine metagenome">
##       <OrganismName>marine metagenome</OrganismName>
##     </Organism>
##   </Description>
##   <Owner>
##     <Name abbreviation="JGI" url="http://jgi.doe.gov/">DOE Joint Genome Institute</Name>
##   </Owner>
##   <Models>
##     <Model>Metagenome or environmental</Model>
##   </Models>
##   <Package display_name="Metagenome or environmental; version 1.0">Metagenome.environmental.1.0</Package>
##   <Attributes>
##     <Attribute attribute_name="collection_date" harmonized_name="collection_date" display_name="collection date">missing</Attribute>
##     <Attribute attribute_name="host" harmonized_name="host" display_name="host">not applicable</Attribute>
##     <Attribute attribute_name="geo_loc_name" harmonized_name="geo_loc_name" display_name="geographic location">missing</Attribute>
##     <Attribute attribute_name="lat_lon" harmonized_name="lat_lon" display_name="latitude and longitude">48.6 N 123.5 W</Attribute>
##     <Attribute attribute_name="isolation_source" harmonized_name="isolation_source" display_name="isolation source">Saanich Inlet, British Columbia, Canada</Attribute>
##   </Attributes>
##   <Links>
##     <Link type="url" label="Gold Stamp ID Gp0091327">https://gold.jgi.doe.gov/projects?id=Gp0091327</Link>
##     <Link type="entrez" target="bioproject" label="PRJNA247822">247822</Link>
##   </Links>
##   <Status status="live" when="2016-06-09T17:54:30.155"/>
## </BioSample>
```



Looks like it's not going to be the easiest time parsing and orienting the data. I found a tutorial online for this....



https://www.stat.berkeley.edu/~statcur/Workshop2/Presentations/XML.pdf




```r
top <- xmlRoot(parsed_biosample_sampledata)
xmlName(top)
```

```
## [1] "para"
```

```r
str(parsed_biosample_sampledata)
```

```
## Classes 'XMLInternalElementNode', 'XMLInternalNode', 'XMLAbstractNode' <externalptr>
```


Find the header 



```r
names(top)
```

```
##   BioSample 
## "BioSample"
```


Give the names of the child nodes of this root



```r
names( top[[1]])
```

```
##           Ids   Description         Owner        Models       Package 
##         "Ids" "Description"       "Owner"      "Models"     "Package" 
##    Attributes         Links        Status 
##  "Attributes"       "Links"      "Status"
```



```r
str( top [[1]])
```

```
## Classes 'XMLInternalElementNode', 'XMLInternalNode', 'XMLAbstractNode' <externalptr>
```



```r
biosample_attributes <- top[[1]] [["Attributes"]]
biosample_attributes
```

```
## <Attributes>
##   <Attribute attribute_name="collection_date" harmonized_name="collection_date" display_name="collection date">missing</Attribute>
##   <Attribute attribute_name="host" harmonized_name="host" display_name="host">not applicable</Attribute>
##   <Attribute attribute_name="geo_loc_name" harmonized_name="geo_loc_name" display_name="geographic location">missing</Attribute>
##   <Attribute attribute_name="lat_lon" harmonized_name="lat_lon" display_name="latitude and longitude">48.6 N 123.5 W</Attribute>
##   <Attribute attribute_name="isolation_source" harmonized_name="isolation_source" display_name="isolation source">Saanich Inlet, British Columbia, Canada</Attribute>
## </Attributes>
```



This is what we are looking for! Except, here is where it is getting a little tricky. I want to be able to specifically take out the latitude and longitude, but it seems that the "Attributes" column is in a weird class of XML that is not as accessible. We may need to get creative as to how we will extract the values. Perhaps a **'grep()'** function?




```r
parse_sampledata <- biosample_attributes["//lat_lon"]
parse_sampledata
```

```
## named list()
## attr(,"class")
## [1] "XMLInternalNodeList" "XMLNodeList"
```


It looks like the 'readHTMLList()' function might be able to help us out here.



```r
biosample_parsed <- readHTMLList(biosample_attributes)
biosample_parsed
```

```
## [1] "missing"                                
## [2] "not applicable"                         
## [3] "missing"                                
## [4] "48.6 N 123.5 W"                         
## [5] "Saanich Inlet, British Columbia, Canada"
```


We got it out! Now lets try to put it into a data frame 

#Put the output into a data frame 



```r
biosample_df <- as.data.frame(biosample_parsed)
biosample_df
```

```
##                          biosample_parsed
## 1                                 missing
## 2                          not applicable
## 3                                 missing
## 4                          48.6 N 123.5 W
## 5 Saanich Inlet, British Columbia, Canada
```




```r
lat_long <- biosample_df[4,]
lat_long
```

```
## [1] 48.6 N 123.5 W
## 4 Levels: 48.6 N 123.5 W missing ... Saanich Inlet, British Columbia, Canada
```


This is the specific data that we are interested in. But is all of the 'Attributes' metadata important?



```r
#Name a rows vector for the dataframe
rows <- c("Collection date", "Host", "Geographic location", "Latitude and longitude", "Isolation source")
rows
```

```
## [1] "Collection date"        "Host"                  
## [3] "Geographic location"    "Latitude and longitude"
## [5] "Isolation source"
```




```r
row.names(biosample_df)<- rows
biosample_df
```

```
##                                               biosample_parsed
## Collection date                                        missing
## Host                                            not applicable
## Geographic location                                    missing
## Latitude and longitude                          48.6 N 123.5 W
## Isolation source       Saanich Inlet, British Columbia, Canada
```

```r
#column.names(biosample_df)<- 
#we should try to make the column names the same as the SRA ID that we are first searching with.  
```



WE'VE GOT IT! Now to clean it up and use more accessions! The next step is to make sure that nothing is 'hard-coded' so that we can reproduce this search with any ID that we put in.
