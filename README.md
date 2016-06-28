# Diana
The DIANA tool for targeted data extraction from DIA mass spectrometry data.
More information at http://quantitativeproteomics.org/diana

## Running Diana
The python scripts and java jar files needed to run Diana can be found in the 'diana' direcory. 

## Compilation
For compilation of the java components of Diana you will need to install maven and Java version 1.8.
The project also depends on components from Anubis, that are currently not in Maven Central. To resolve these dependencies, first clone the Anubis repository (http://github.com/fickludd/anubis) and run the python installation script provided with Anubis. Then compile Diana.
Diana components are compiled using maven by decending into the respective directories and running 'mvn install'. This should be done in the following order:
DianaExtractor
DianaAnalyzeLib
DianaScorer 

### License
The Diana code is licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
