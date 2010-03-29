package org.fhcrc.cpl.toolbox.chem;

import java.util.HashMap;

/**
 * Created by IntelliJ IDEA.
 * User: tholzman
 * Date: Mar 25, 2010
 * Time: 12:08:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class Elements {
    public static HashMap<String, Element> Elements = new HashMap<String, Element>();

    static {
        Elements.put("H", new Element("H",1.007825032, 1));
        Elements.put("He", new Element("He",4.002603254, 2));
        Elements.put("Li", new Element("Li",7.01600455, 3));
        Elements.put("Be", new Element("Be",9.0121822, 4));
        Elements.put("B", new Element("B",11.0093054, 5));
        Elements.put("C", new Element("C",12D, 6));
        Elements.put("N", new Element("N",14.003074, 7));
        Elements.put("O", new Element("O",15.99491462, 8));
        Elements.put("F", new Element("F",18.99840322, 9));
        Elements.put("Ne", new Element("Ne",19.99244018, 10));
        Elements.put("Na", new Element("Na",22.98976928, 11));
        Elements.put("Mg", new Element("Mg",23.9850417, 12));
        Elements.put("Al", new Element("Al",26.98153863, 13));
        Elements.put("Si", new Element("Si",27.97692653, 14));
        Elements.put("P", new Element("P",30.97376163, 15));
        Elements.put("S", new Element("S",31.972071, 16));
        Elements.put("Cl", new Element("Cl",34.96885268, 17));
        Elements.put("Ar", new Element("Ar",39.96238312, 18));
        Elements.put("K", new Element("K",38.96370668, 19));
        Elements.put("Ca", new Element("Ca",39.96259098, 20));
        Elements.put("Sc", new Element("Sc",44.9559119, 21));
        Elements.put("Ti", new Element("Ti",47.9479463, 22));
        Elements.put("V", new Element("V",50.9439595, 23));
        Elements.put("Cr", new Element("Cr",51.9405075, 24));
        Elements.put("Mn", new Element("Mn",54.9380451, 25));
        Elements.put("Fe", new Element("Fe",55.9349375, 26));
        Elements.put("Co", new Element("Co",58.933195, 27));
        Elements.put("Ni", new Element("Ni",57.9353429, 28));
        Elements.put("Cu", new Element("Cu",62.9295975, 29));
        Elements.put("Zn", new Element("Zn",63.9291422, 30));
        Elements.put("Ga", new Element("Ga",68.9255736, 31));
        Elements.put("Ge", new Element("Ge",73.9211778, 32));
        Elements.put("As", new Element("As",74.9215965, 33));
        Elements.put("Se", new Element("Se",79.9165213, 34));
        Elements.put("Br", new Element("Br",78.9183371, 35));
        Elements.put("Kr", new Element("Kr",83.911507, 36));
        Elements.put("Rb", new Element("Rb",84.91178974, 37));
        Elements.put("Sr", new Element("Sr",87.9056121, 38));
        Elements.put("Y", new Element("Y",88.9058483, 39));
        Elements.put("Zr", new Element("Zr",89.9047044, 40));
        Elements.put("Nb", new Element("Nb",92.9063781, 41));
        Elements.put("Mo", new Element("Mo",97.9054082, 42));
        Elements.put("Tc", new Element("Tc", 98D, 43));
        Elements.put("Ru", new Element("Ru",101.9043493, 44));
        Elements.put("Rh", new Element("Rh",102.905504, 45));
        Elements.put("Pd", new Element("Pd",105.903486, 46));
        Elements.put("Ag", new Element("Ag",106.905097, 47));
        Elements.put("Cd", new Element("Cd",113.9033585, 48));
        Elements.put("In", new Element("In",114.903878, 49));
        Elements.put("Sn", new Element("Sn",119.9021947, 50));
        Elements.put("Sb", new Element("Sb",120.9038157, 51));
        Elements.put("Te", new Element("Te",129.9062244, 52));
        Elements.put("I", new Element("I",126.904473, 53));
        Elements.put("Xe", new Element("Xe",131.9041535, 54));
        Elements.put("Cs", new Element("Cs",132.9054519, 55));
        Elements.put("Ba", new Element("Ba",137.9052472, 56));
        Elements.put("La", new Element("La",138.9063533, 57));
        Elements.put("Ce", new Element("Ce",139.9054387, 58));
        Elements.put("Pr", new Element("Pr",140.9076528, 59));
        Elements.put("Nd", new Element("Nd",141.9077233, 60));
        Elements.put("Pm", new Element("Pm", 145D, 61));
        Elements.put("Sm", new Element("Sm",151.9197324, 62));
        Elements.put("Eu", new Element("Eu",152.9212303, 63));
        Elements.put("Gd", new Element("Gd",157.9241039, 64));
        Elements.put("Tb", new Element("Tb",158.9253468, 65));
        Elements.put("Dy", new Element("Dy",163.9291748, 66));
        Elements.put("Ho", new Element("Ho",164.9303221, 67));
        Elements.put("Er", new Element("Er",165.9302931, 68));
        Elements.put("Tm", new Element("Tm",168.9342133, 69));
        Elements.put("Yb", new Element("Yb",173.9388621, 70));
        Elements.put("Lu", new Element("Lu",174.9407718, 71));
        Elements.put("Hf", new Element("Hf",179.94655, 72));
        Elements.put("Ta", new Element("Ta",180.9479958, 73));
        Elements.put("W", new Element("W",183.9509312, 74));
        Elements.put("Re", new Element("Re",186.9557531, 75));
        Elements.put("Os", new Element("Os",191.9614807, 76));
        Elements.put("Ir", new Element("Ir",192.9629264, 77));
        Elements.put("Pt", new Element("Pt",194.9647911, 78));
        Elements.put("Au", new Element("Au",196.9665687, 79));
        Elements.put("Hg", new Element("Hg",201.970643, 80));
        Elements.put("Tl", new Element("Tl",204.9744275, 81));
        Elements.put("Pb", new Element("Pb",207.9766521, 82));
        Elements.put("Bi", new Element("Bi",208.9803987, 83));
        Elements.put("Po", new Element("Po", 209D, 84));
        Elements.put("At", new Element("At", 210D, 85));
        Elements.put("Rn", new Element("Rn", 222D, 86));
        Elements.put("Fr", new Element("Fr", 223D, 87));
        Elements.put("Ra", new Element("Ra",226.0254, 88));
        Elements.put("Ac", new Element("Ac", 227D, 89));
        Elements.put("Th", new Element("Th",232.0380553, 90));
        Elements.put("Pa", new Element("Pa",231.035884, 91));
        Elements.put("U", new Element("U",238.0507882, 92));
        Elements.put("Np", new Element("Np",237.0482, 93));
        Elements.put("Pu", new Element("Pu", 244D, 94));
        Elements.put("Am", new Element("Am", 243D, 95));
        Elements.put("Cm", new Element("Cm", 247D, 96));
        Elements.put("Bk", new Element("Bk", 247D, 97));
        Elements.put("Cf", new Element("Cf", 251D, 98));
        Elements.put("Es", new Element("Es", 252D, 99));
        Elements.put("Fm", new Element("Fm", 257D, 100));
        Elements.put("Md", new Element("Md", 258D, 101));
        Elements.put("No", new Element("No", 259D, 102));
        Elements.put("Lr", new Element("Lr", 260D, 103));
        Elements.put("Rf", new Element("Rf", 261D, 104));
        Elements.put("Db", new Element("Db", 262D, 105));
        Elements.put("Sg", new Element("Sg", 263D, 106));
        Elements.put("Bh", new Element("Bh", 262D, 107));
        Elements.put("Hs", new Element("Hs", 265D, 108));
        Elements.put("Mt", new Element("Mt", 266D, 109));

    }
}
