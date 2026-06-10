CREATE DATABASE mitoCancer;

CREATE TABLE mitoCartaMaster (
geneId SERIAL PRIMARY KEY,
entrezId INT UNIQUE,
geneSymbol VARCHAR(50) NOT NULL,
mouseOrtholog VARCHAR(50),
geneDescription TEXT,
uniprotId VARCHAR(25),
mitoPathway TEXT,
mitocartaScore FLOAT NOT NULL,
targetpScore FLOAT,
hg19chr VARCHAR(10)
inMitocarta BOOLEAN DEFAULT TRUE,
);

CREATE TABLE ensemblIdMapping(
ensemblgeneId VARCHAR(50) PRIMARY KEY,
geneId INT REFERENCES mitoCartaMaster(geneId),
idRank INT NOT NULL DEFAULT 1
);