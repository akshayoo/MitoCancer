CREATE DATABASE mitoCancer;

CREATE TABLE mitoCartaMaster (
geneId SERIAL PRIMARY KEY,
entrezId INT UNIQUE,
geneSymbol VARCHAR(50) NOT NULL,
mouseOrtholog VARCHAR(50),
geneDescription TEXT,
submitoLocalization VARCHAR(25),
uniprotId VARCHAR(25),
mitoPathway TEXT,
mitocartaScore FLOAT,
targetpScore FLOAT,
hg19chr VARCHAR(10),
inMitocarta BOOLEAN DEFAULT TRUE
);

CREATE TABLE ensemblIdMapping(
ensemblgeneId VARCHAR(50) PRIMARY KEY,
geneId INT REFERENCES mitoCartaMaster(geneId),
idRank INT NOT NULL DEFAULT 1
);

CREATE TABLE uniprotMitoData (
uniprotId VARCHAR(20) PRIMARY KEY,
geneId INT REFERENCES mitocartamaster(geneId),
proteinEntry VARCHAR(25),
proteinName TEXT,
proteinExistence VARCHAR(50),
subcellularLocation TEXT,
sequenceLength INT,
molWeight INT,
isCanonical BOOLEAN DEFAULT TRUE,
dateRetrieved TIMESTAMP DEFAULT NOW()
);

CREATE TABLE uniprotMitoDiseaseData (
id SERIAL PRIMARY KEY,
uniprotId VARCHAR(20) REFERENCES uniprotmitodata(uniprotid),
geneId INT REFERENCES mitocartamaster(geneId),
diseaseId VARCHAR(100),
diseaseAcronym VARCHAR(20), 
mimId VARCHAR(10), 
diseaseDescription TEXT,    
isCancer BOOLEAN       
);