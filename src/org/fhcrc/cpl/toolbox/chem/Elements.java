/*
 * Copyright (c) 2003-2012 Fred Hutchinson Cancer Research Center
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.fhcrc.cpl.toolbox.chem;

import org.apache.log4j.Logger;
import org.fhcrc.cpl.toolbox.proteomics.feature.filehandler.APMLFeatureFileHandler;

import java.util.*;

/**
 * Provides information about all of the elements, by instantiating Element objects.
 *
 * In particular, instantiates a HashMap of Elements, mapped by symbol.  Each Element has
 * its monoisotopic and average masses defined.  For elements that we "care about",
 * peak masses and frequencies are defined.  Those are left null for elements we don't "care about"
 *
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 25, 2010
 * Time: 12:08:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class Elements
{
    protected static Logger _log = Logger.getLogger(APMLFeatureFileHandler.class);
    protected static List<String> elementsToLog = new ArrayList<String>(
            Arrays.asList( new String[] { "C","H","O","N","P","S","Cl" }));

    protected static HashMap<String, Element> ELEMENTS_BY_SYMBOL = new HashMap<String, Element>();

    //todo: element should store most common atomic weight
    static
    {
        ELEMENTS_BY_SYMBOL.put("H",  new Element("H",1,new double[] {1.007825032,2.014101778}, new double[] {0.999885,1.15E-4}));
        ELEMENTS_BY_SYMBOL.put("He",  new Element("He",2,new double[] {3.016029319,4.002603254}, new double[] {1.0E-6,0.999999}));
        ELEMENTS_BY_SYMBOL.put("Li",  new Element("Li",3,new double[] {6.015122795,7.01600455}, new double[] {0.0759,0.9241}));
        ELEMENTS_BY_SYMBOL.put("Be",  new Element("Be",4,new double[] {9.0121822}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("B",  new Element("B",5,new double[] {10.012937,11.0093054}, new double[] {0.199,0.801}));
        ELEMENTS_BY_SYMBOL.put("C",  new Element("C",6,new double[] {12.0,13.00335484}, new double[] {0.9893,0.0107}));
        ELEMENTS_BY_SYMBOL.put("N",  new Element("N",7,new double[] {14.003074,15.0001089}, new double[] {0.99632,0.00368}));
        ELEMENTS_BY_SYMBOL.put("O",  new Element("O",8,new double[] {15.99491462,16.9991317,17.999161}, new double[] {0.99757,3.8E-4,0.00205}));
        ELEMENTS_BY_SYMBOL.put("F",  new Element("F",9,new double[] {18.99840322}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Ne",  new Element("Ne",10,new double[] {19.99244018,20.99384668,21.99138511}, new double[] {0.9048,0.0027,0.0925}));
        ELEMENTS_BY_SYMBOL.put("Na",  new Element("Na",11,new double[] {22.98976928}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Mg",  new Element("Mg",12,new double[] {23.9850417,24.98583692,25.98259293}, new double[] {0.7899,0.1,0.1101}));
        ELEMENTS_BY_SYMBOL.put("Al",  new Element("Al",13,new double[] {26.98153863}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Si",  new Element("Si",14,new double[] {27.97692653,28.9764947,29.97377017}, new double[] {0.922297,0.046832,0.030872}));
        ELEMENTS_BY_SYMBOL.put("P",  new Element("P",15,new double[] {30.97376163}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("S",  new Element("S",16,new double[] {31.972071,32.97145876,33.9678669,35.96708076}, new double[] {0.9493,0.0076,0.0429,2.0E-4}));
        ELEMENTS_BY_SYMBOL.put("Cl",  new Element("Cl",17,new double[] {34.96885268,36.96590259}, new double[] {0.7578,0.2422}));
        ELEMENTS_BY_SYMBOL.put("Ar",  new Element("Ar",18,new double[] {35.96754511,37.9627324,39.96238312}, new double[] {0.003365,6.32E-4,0.996003}));
        ELEMENTS_BY_SYMBOL.put("K",  new Element("K",19,new double[] {38.96370668,39.96399848,40.96182576}, new double[] {0.932581,1.17E-4,0.067302}));
        ELEMENTS_BY_SYMBOL.put("Ca",  new Element("Ca",20,new double[] {39.96259098,41.95861801,42.9587666,43.9554818,45.9536926,47.952534}, new double[] {0.96941,0.00647,0.00135,0.02086,4.0E-5,0.00187}));
        ELEMENTS_BY_SYMBOL.put("Sc",  new Element("Sc",21,new double[] {44.9559119}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Ti",  new Element("Ti",22,new double[] {45.9526316,46.9517631,47.9479463,48.94787,49.9447912}, new double[] {0.0825,0.0744,0.7372,0.0541,0.0518}));
        ELEMENTS_BY_SYMBOL.put("V",  new Element("V",23,new double[] {49.9471585,50.9439595}, new double[] {0.0025,0.9975}));
        ELEMENTS_BY_SYMBOL.put("Cr",  new Element("Cr",24,new double[] {49.9460442,51.9405075,52.9406494,53.9388804}, new double[] {0.04345,0.83789,0.09501,0.02365}));
        ELEMENTS_BY_SYMBOL.put("Mn",  new Element("Mn",25,new double[] {54.9380451}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Fe",  new Element("Fe",26,new double[] {53.9396105,55.9349375,56.935394,57.9332756}, new double[] {0.05845,0.91754,0.02119,0.00282}));
        ELEMENTS_BY_SYMBOL.put("Co",  new Element("Co",27,new double[] {58.933195}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Ni",  new Element("Ni",28,new double[] {57.9353429,59.9307864,60.931056,61.9283451,63.927966}, new double[] {0.680769,0.262231,0.011399,0.036345,0.009256}));
        ELEMENTS_BY_SYMBOL.put("Cu",  new Element("Cu",29,new double[] {62.9295975,64.9277895}, new double[] {0.6917,0.3083}));
        ELEMENTS_BY_SYMBOL.put("Zn",  new Element("Zn",30,new double[] {63.9291422,65.9260334,66.9271273,67.9248442,69.9253193}, new double[] {0.4863,0.279,0.041,0.1875,0.0062}));
        ELEMENTS_BY_SYMBOL.put("Ga",  new Element("Ga",31,new double[] {68.9255736,70.9247013}, new double[] {0.60108,0.39892}));
        ELEMENTS_BY_SYMBOL.put("Ge",  new Element("Ge",32,new double[] {69.9242474,71.9220758,72.9234589,73.9211778,75.9214026}, new double[] {0.2084,0.2754,0.0773,0.3628,0.0761}));
        ELEMENTS_BY_SYMBOL.put("As",  new Element("As",33,new double[] {74.9215965}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Se",  new Element("Se",34,new double[] {73.9224764,75.9192136,76.919914,77.9173091,79.9165213,81.9166994}, new double[] {0.0089,0.0937,0.0763,0.2377,0.4961,0.0873}));
        ELEMENTS_BY_SYMBOL.put("Br",  new Element("Br",35,new double[] {78.9183371,80.9162906}, new double[] {0.5069,0.4931}));
        ELEMENTS_BY_SYMBOL.put("Kr",  new Element("Kr",36,new double[] {77.9203648,79.916379,81.9134836,82.914136,83.911507,85.91061073}, new double[] {0.0035,0.0228,0.1158,0.1149,0.57,0.173}));
        ELEMENTS_BY_SYMBOL.put("Rb",  new Element("Rb",37,new double[] {84.91178974,86.90918053}, new double[] {0.7217,0.2783}));
        ELEMENTS_BY_SYMBOL.put("Sr",  new Element("Sr",38,new double[] {83.913425,85.9092602,86.9088771,87.9056121}, new double[] {0.0056,0.0986,0.07,0.8258}));
        ELEMENTS_BY_SYMBOL.put("Y",  new Element("Y",39,new double[] {88.9058483}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Zr",  new Element("Zr",40,new double[] {89.9047044,90.9056458,91.9050408,93.9063152,95.9082734}, new double[] {0.5145,0.1122,0.1715,0.1738,0.028}));
        ELEMENTS_BY_SYMBOL.put("Nb",  new Element("Nb",41,new double[] {92.9063781}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Mo",  new Element("Mo",42,new double[] {91.906811,93.9050883,94.9058421,95.9046795,96.9060215,97.9054082,99.907477}, new double[] {0.1484,0.0925,0.1592,0.1668,0.0955,0.2413,0.0963}));
        ELEMENTS_BY_SYMBOL.put("Ru",  new Element("Ru",44,new double[] {95.907598,97.905287,98.9059393,99.9042195,100.9055821,101.9043493,103.905433}, new double[] {0.0554,0.0187,0.1276,0.126,0.1706,0.3155,0.1862}));
        ELEMENTS_BY_SYMBOL.put("Rh",  new Element("Rh",45,new double[] {102.905504}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Pd",  new Element("Pd",46,new double[] {101.905609,103.904036,104.905085,105.903486,107.903892,109.905153}, new double[] {0.0102,0.1114,0.2233,0.2733,0.2646,0.1172}));
        ELEMENTS_BY_SYMBOL.put("Ag",  new Element("Ag",47,new double[] {106.905097,108.904752}, new double[] {0.51839,0.48161}));
        ELEMENTS_BY_SYMBOL.put("Cd",  new Element("Cd",48,new double[] {105.906459,107.904184,109.9030021,110.9041781,111.9027578,112.9044017,113.9033585,115.904756}, new double[] {0.0125,0.0089,0.1249,0.128,0.2413,0.1222,0.2873,0.0749}));
        ELEMENTS_BY_SYMBOL.put("In",  new Element("In",49,new double[] {112.904058,114.903878}, new double[] {0.0429,0.9571}));
        ELEMENTS_BY_SYMBOL.put("Sn",  new Element("Sn",50,new double[] {111.904818,113.902779,114.903342,115.901741,116.902952,117.901603,118.903308,119.9021947,121.903439,123.9052739}, new double[] {0.0097,0.0066,0.0034,0.1454,0.0768,0.2422,0.0859,0.3258,0.0463,0.0579}));
        ELEMENTS_BY_SYMBOL.put("Sb",  new Element("Sb",51,new double[] {120.9038157,122.904214}, new double[] {0.5721,0.4279}));
        ELEMENTS_BY_SYMBOL.put("Te",  new Element("Te",52,new double[] {119.90402,121.9030439,122.90427,123.9028179,124.9044307,125.9033117,127.9044631,129.9062244}, new double[] {9.0E-4,0.0255,0.0089,0.0474,0.0707,0.1884,0.3174,0.3408}));
        ELEMENTS_BY_SYMBOL.put("I",  new Element("I",53,new double[] {126.904473}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Xe",  new Element("Xe",54,new double[] {123.905893,125.904274,127.9035313,128.9047794,129.903508,130.9050824,131.9041535,133.9053945,135.907219}, new double[] {9.0E-4,9.0E-4,0.0192,0.2644,0.0408,0.2118,0.2689,0.1044,0.0887}));
        ELEMENTS_BY_SYMBOL.put("Cs",  new Element("Cs",55,new double[] {132.9054519}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Ba",  new Element("Ba",56,new double[] {129.9063208,131.9050613,133.9045084,134.9056886,135.9045759,136.9058274,137.9052472}, new double[] {0.00106,0.00101,0.02417,0.06592,0.07854,0.11232,0.71698}));
        ELEMENTS_BY_SYMBOL.put("La",  new Element("La",57,new double[] {137.907112,138.9063533}, new double[] {9.0E-4,0.9991}));
        ELEMENTS_BY_SYMBOL.put("Ce",  new Element("Ce",58,new double[] {135.907172,137.905991,139.9054387,141.909244}, new double[] {0.00185,0.00251,0.8845,0.11114}));
        ELEMENTS_BY_SYMBOL.put("Pr",  new Element("Pr",59,new double[] {140.9076528}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Nd",  new Element("Nd",60,new double[] {141.9077233,142.9098143,143.9100873,144.9125736,145.9131169,147.916893,149.920891}, new double[] {0.272,0.122,0.238,0.083,0.172,0.057,0.056}));
        ELEMENTS_BY_SYMBOL.put("Sm",  new Element("Sm",62,new double[] {143.911999,146.9148979,147.9148227,148.9171847,149.9172755,151.9197324,153.9222093}, new double[] {0.0307,0.1499,0.1124,0.1382,0.0738,0.2675,0.2275}));
        ELEMENTS_BY_SYMBOL.put("Eu",  new Element("Eu",63,new double[] {150.9198502,152.9212303}, new double[] {0.4781,0.5219}));
        ELEMENTS_BY_SYMBOL.put("Gd",  new Element("Gd",64,new double[] {151.919791,153.9208656,154.922622,155.9221227,156.9239601,157.9241039,159.9270541}, new double[] {0.0020,0.0218,0.148,0.2047,0.1565,0.2484,0.2186}));
        ELEMENTS_BY_SYMBOL.put("Tb",  new Element("Tb",65,new double[] {158.9253468}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Dy",  new Element("Dy",66,new double[] {155.924283,157.924409,159.9251975,160.9269334,161.9267984,162.9287312,163.9291748}, new double[] {6.0E-4,0.0010,0.0234,0.1891,0.2551,0.249,0.2818}));
        ELEMENTS_BY_SYMBOL.put("Ho",  new Element("Ho",67,new double[] {164.9303221}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Er",  new Element("Er",68,new double[] {161.928778,163.9292,165.9302931,166.9320482,167.9323702,169.9354643}, new double[] {0.0014,0.0161,0.3361,0.2293,0.2678,0.1493}));
        ELEMENTS_BY_SYMBOL.put("Tm",  new Element("Tm",69,new double[] {168.9342133}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Yb",  new Element("Yb",70,new double[] {167.933897,169.9347618,170.9363258,171.9363815,172.9382108,173.9388621,175.9425717}, new double[] {0.0013,0.0304,0.1428,0.2183,0.1613,0.3183,0.1276}));
        ELEMENTS_BY_SYMBOL.put("Lu",  new Element("Lu",71,new double[] {174.9407718,175.9426863}, new double[] {0.9741,0.0259}));
        ELEMENTS_BY_SYMBOL.put("Hf",  new Element("Hf",72,new double[] {173.940046,175.9414086,176.9432207,177.9436988,178.9458161,179.94655}, new double[] {0.0016,0.0526,0.186,0.2728,0.1362,0.3508}));
        ELEMENTS_BY_SYMBOL.put("Ta",  new Element("Ta",73,new double[] {179.9474648,180.9479958}, new double[] {1.2E-4,0.99988}));
        ELEMENTS_BY_SYMBOL.put("W",  new Element("W",74,new double[] {179.946704,181.9482042,182.950223,183.9509312,185.9543641}, new double[] {0.0012,0.265,0.1431,0.3064,0.2843}));
        ELEMENTS_BY_SYMBOL.put("Re",  new Element("Re",75,new double[] {184.952955,186.9557531}, new double[] {0.374,0.626}));
        ELEMENTS_BY_SYMBOL.put("Os",  new Element("Os",76,new double[] {183.9524891,185.9538382,186.9557505,187.9558382,188.9581475,189.958447,191.9614807}, new double[] {2.0E-4,0.0159,0.0196,0.1324,0.1615,0.2626,0.4078}));
        ELEMENTS_BY_SYMBOL.put("Ir",  new Element("Ir",77,new double[] {190.960594,192.9629264}, new double[] {0.373,0.627}));
        ELEMENTS_BY_SYMBOL.put("Pt",  new Element("Pt",78,new double[] {189.959932,191.961038,193.9626803,194.9647911,195.9649515,197.967893}, new double[] {1.4E-4,0.00782,0.32967,0.33832,0.25242,0.07163}));
        ELEMENTS_BY_SYMBOL.put("Au",  new Element("Au",79,new double[] {196.9665687}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Hg",  new Element("Hg",80,new double[] {195.965833,197.966769,198.9682799,199.968326,200.9703023,201.970643,203.9734939}, new double[] {0.0015,0.0997,0.1687,0.231,0.1318,0.2986,0.0687}));
        ELEMENTS_BY_SYMBOL.put("Tl",  new Element("Tl",81,new double[] {202.9723442,204.9744275}, new double[] {0.29524,0.70476}));
        ELEMENTS_BY_SYMBOL.put("Pb",  new Element("Pb",82,new double[] {203.9730436,205.9744653,206.9758969,207.9766521}, new double[] {0.014,0.241,0.221,0.524}));
        ELEMENTS_BY_SYMBOL.put("Bi",  new Element("Bi",83,new double[] {208.9803987}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Th",  new Element("Th",90,new double[] {232.0380553}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("Pa",  new Element("Pa",91,new double[] {231.035884}, new double[] {1.0}));
        ELEMENTS_BY_SYMBOL.put("U",  new Element("U",92,new double[] {234.0409521,235.0439299,238.0507882}, new double[] {5.5E-5,0.0072,0.992745}));


        if (_log.isDebugEnabled())
        {
            for (String element : elementsToLog)
                _log.debug(get(element));
        }



//        ELEMENTS.put("Db", new Element("Db",262.0,262.0,105));
//        ELEMENTS.put("Cs", new Element("Cs",132.9054519,132.90543,55));
//        ELEMENTS.put("Cu", new Element("Cu",62.9295975,63.546,29));
//        ELEMENTS.put("Kr", new Element("Kr",83.911507,83.8,36));
//        ELEMENTS.put("Cl", new Element("Cl",34.96885268,35.4527,17));
//        ELEMENTS.put("Cm", new Element("Cm",247.0,247.0,96));
//        ELEMENTS.put("Co", new Element("Co",58.933195,58.9332,27));
//        ELEMENTS.put("Cr", new Element("Cr",51.9405075,51.9961,24));
//        ELEMENTS.put("Li", new Element("Li",7.01600455,6.941,3));
//        ELEMENTS.put("Cd", new Element("Cd",113.9033585,112.411,48));
//        ELEMENTS.put("Cf", new Element("Cf",251.0,251.0,98));
//        ELEMENTS.put("Ce", new Element("Ce",139.9054387,140.115,58));
//        ELEMENTS.put("La", new Element("La",138.9063533,138.9055,57));
//        ELEMENTS.put("Lu", new Element("Lu",174.9407718,174.967,71));
//        ELEMENTS.put("Tl", new Element("Tl",204.9744275,204.3833,81));
//        ELEMENTS.put("Tm", new Element("Tm",168.9342133,168.93421,69));
//        ELEMENTS.put("Th", new Element("Th",232.0380553,232.0381,90));
//        ELEMENTS.put("Ti", new Element("Ti",47.9479463,47.88,22));
//        ELEMENTS.put("Te", new Element("Te",129.9062244,127.6,52));
//        ELEMENTS.put("Lr", new Element("Lr",260.0,260.0,103));
//        ELEMENTS.put("Dy", new Element("Dy",163.9291748,162.5,66));
//        ELEMENTS.put("Ta", new Element("Ta",180.9479958,180.9479,73));
//        ELEMENTS.put("Tc", new Element("Tc",98.0,98.0,43));
//        ELEMENTS.put("Mg", new Element("Mg",23.9850417,24.305,12));
//        ELEMENTS.put("Tb", new Element("Tb",158.9253468,158.92534,65));
//        ELEMENTS.put("Md", new Element("Md",258.0,258.0,101));
//        ELEMENTS.put("F", new Element("F",18.99840322,18.9984032,9));
//        ELEMENTS.put("Fe", new Element("Fe",55.9349375,55.847,26));
//        ELEMENTS.put("B", new Element("B",11.0093054,10.811,5));
//        ELEMENTS.put("C", new Element("C",12.0,12.011,6));
//        ELEMENTS.put("Mt", new Element("Mt",266.0,266.0,109));
//        ELEMENTS.put("N", new Element("N",14.003074,14.00674,7));
//        ELEMENTS.put("O", new Element("O",15.99491462,15.9994,8));
//        ELEMENTS.put("H", new Element("H",1.007825032,1.00794,1));
//        ELEMENTS.put("Eu", new Element("Eu",152.9212303,151.965,63));
//        ELEMENTS.put("Mo", new Element("Mo",97.9054082,95.94,42));
//        ELEMENTS.put("I", new Element("I",126.904473,126.90447,53));
//        ELEMENTS.put("Mn", new Element("Mn",54.9380451,54.93805,25));
//        ELEMENTS.put("K", new Element("K",38.96370668,39.0983,19));
//        ELEMENTS.put("Er", new Element("Er",165.9302931,167.26,68));
//        ELEMENTS.put("U", new Element("U",238.0507882,238.0289,92));
//        ELEMENTS.put("W", new Element("W",183.9509312,183.85,74));
//        ELEMENTS.put("Es", new Element("Es",252.0,252.0,99));
//        ELEMENTS.put("V", new Element("V",50.9439595,50.9415,23));
//        ELEMENTS.put("Ni", new Element("Ni",57.9353429,58.6934,28));
//        ELEMENTS.put("P", new Element("P",30.97376163,30.973762,15));
//        ELEMENTS.put("S", new Element("S",31.972071,32.066,16));
//        ELEMENTS.put("Nd", new Element("Nd",141.9077233,144.24,60));
//        ELEMENTS.put("Ne", new Element("Ne",19.99244018,20.1797,10));
//        ELEMENTS.put("Nb", new Element("Nb",92.9063781,92.90638,41));
//        ELEMENTS.put("Y", new Element("Y",88.9058483,88.90585,39));
//        ELEMENTS.put("Na", new Element("Na",22.98976928,22.989768,11));
//        ELEMENTS.put("Ge", new Element("Ge",73.9211778,72.61,32));
//        ELEMENTS.put("Gd", new Element("Gd",157.9241039,157.25,64));
//        ELEMENTS.put("Ga", new Element("Ga",68.9255736,69.723,31));
//        ELEMENTS.put("No", new Element("No",259.0,259.0,102));
//        ELEMENTS.put("Np", new Element("Np",237.0482,237.0482,93));
//        ELEMENTS.put("Fr", new Element("Fr",223.0,223.0,87));
//        ELEMENTS.put("Fm", new Element("Fm",257.0,257.0,100));
//        ELEMENTS.put("Yb", new Element("Yb",173.9388621,173.04,70));
//        ELEMENTS.put("Pt", new Element("Pt",194.9647911,195.08,78));
//        ELEMENTS.put("Pu", new Element("Pu",244.0,244.0,94));
//        ELEMENTS.put("Pr", new Element("Pr",140.9076528,140.90765,59));
//        ELEMENTS.put("Hg", new Element("Hg",201.970643,200.59,80));
//        ELEMENTS.put("Hf", new Element("Hf",179.94655,178.49,72));
//        ELEMENTS.put("He", new Element("He",4.002603254,4.002602,2));
//        ELEMENTS.put("Pd", new Element("Pd",105.903486,106.42,46));
//        ELEMENTS.put("Pa", new Element("Pa",231.035884,213.0359,91));
//        ELEMENTS.put("Ho", new Element("Ho",164.9303221,164.93032,67));
//        ELEMENTS.put("Pb", new Element("Pb",207.9766521,207.2,82));
//        ELEMENTS.put("Pm", new Element("Pm",145.0,145.0,61));
//        ELEMENTS.put("Po", new Element("Po",209.0,209.0,84));
//        ELEMENTS.put("Hs", new Element("Hs",265.0,265.0,108));
//        ELEMENTS.put("Xe", new Element("Xe",131.9041535,131.29,54));
//        ELEMENTS.put("Os", new Element("Os",191.9614807,190.2,76));
//        ELEMENTS.put("Se", new Element("Se",79.9165213,78.96,34));
//        ELEMENTS.put("Au", new Element("Au",196.9665687,196.96654,79));
//        ELEMENTS.put("In", new Element("In",114.903878,114.82,49));
//        ELEMENTS.put("Sc", new Element("Sc",44.9559119,44.95591,21));
//        ELEMENTS.put("Ar", new Element("Ar",39.96238312,39.948,18));
//        ELEMENTS.put("Si", new Element("Si",27.97692653,28.0855,14));
//        ELEMENTS.put("At", new Element("At",210.0,210.0,85));
//        ELEMENTS.put("Sg", new Element("Sg",263.0,263.0,106));
//        ELEMENTS.put("As", new Element("As",74.9215965,74.92159,33));
//        ELEMENTS.put("Sn", new Element("Sn",119.9021947,118.71,50));
//        ELEMENTS.put("Sm", new Element("Sm",151.9197324,150.36,62));
//        ELEMENTS.put("Ba", new Element("Ba",137.9052472,137.327,56));
//        ELEMENTS.put("Sr", new Element("Sr",87.9056121,87.62,38));
//        ELEMENTS.put("Ir", new Element("Ir",192.9629264,192.22,77));
//        ELEMENTS.put("Ru", new Element("Ru",101.9043493,101.07,44));
//        ELEMENTS.put("Ag", new Element("Ag",106.905097,107.8682,47));
//        ELEMENTS.put("Ac", new Element("Ac",227.0,227.0,89));
//        ELEMENTS.put("Am", new Element("Am",243.0,243.0,95));
//        ELEMENTS.put("Sb", new Element("Sb",120.9038157,121.757,51));
//        ELEMENTS.put("Al", new Element("Al",26.98153863,26.981539,13));
//        ELEMENTS.put("Rb", new Element("Rb",84.91178974,85.4678,37));
//        ELEMENTS.put("Re", new Element("Re",186.9557531,186.207,75));
//        ELEMENTS.put("Rf", new Element("Rf",261.0,261.0,104));
//        ELEMENTS.put("Rh", new Element("Rh",102.905504,102.9055,45));
//        ELEMENTS.put("Br", new Element("Br",78.9183371,79.904,35));
//        ELEMENTS.put("Ca", new Element("Ca",39.96259098,40.078,20));
//        ELEMENTS.put("Rn", new Element("Rn",222.0,222.0,86));
//        ELEMENTS.put("Bh", new Element("Bh",262.0,262.0,107));
//        ELEMENTS.put("Bi", new Element("Bi",208.9803987,208.98037,83));
//        ELEMENTS.put("Be", new Element("Be",9.0121822,9.012182,4));
//        ELEMENTS.put("Zn", new Element("Zn",63.9291422,65.39,30));
//        ELEMENTS.put("Zr", new Element("Zr",89.9047044,91.224,40));
//        ELEMENTS.put("Ra", new Element("Ra",226.0254,226.0254,88));
//        ELEMENTS.put("Bk", new Element("Bk",247.0,247.0,97));
//
//
//        //define the isotopic peak masses and frequencies for the elements that we care about
//        ELEMENTS.get("C").setIsotopicPeakMasses(new double[] { 12.0, 13.00335 });
//        ELEMENTS.get("C").setIsotopicPeakFrequencies(new double[] { .9893, .0107 });
//
//        ELEMENTS.get("H").setIsotopicPeakMasses(new double[] { 1.007825, 2.014102 });
//        ELEMENTS.get("H").setIsotopicPeakFrequencies(new double[] { .999885, .000115 });
//
//        ELEMENTS.get("N").setIsotopicPeakMasses(new double[] { 14.003074, 15.000109 });
//        ELEMENTS.get("N").setIsotopicPeakFrequencies(new double[] { .99632, .00368 });
//
//        ELEMENTS.get("O").setIsotopicPeakMasses(new double[] { 15.994915, 16.999131, 17.999159 });
//        ELEMENTS.get("O").setIsotopicPeakFrequencies(new double[] { .99757, 0.00038, .00205 });
//
//        ELEMENTS.get("S").setIsotopicPeakMasses(new double[] { 31.972072, 32.971459, 33.967868, 34.9674, 35.967079 });
//        ELEMENTS.get("S").setIsotopicPeakFrequencies(new double[] { .9493, 0.0076, .0420, 0, .00011 });
//
//        ELEMENTS.get("P").setIsotopicPeakMasses(new double[] { 30.973762 });
//        ELEMENTS.get("P").setIsotopicPeakFrequencies(new double[] { 1 });
//
//        ELEMENTS.get("Cl").setIsotopicPeakMasses(new double[] { 34.96885272, 35.96830698, 36.96590259 } );
//        ELEMENTS.get("Cl").setIsotopicPeakFrequencies(new double[] { .7578, 0, .2422 } );
    }

    public static Element get(String elementName)
    {
        return ELEMENTS_BY_SYMBOL.get(elementName);
    }
}
