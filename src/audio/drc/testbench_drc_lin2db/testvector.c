// SPDX-License-Identifier: BSD-3-Clause
//
// Copyright(c) 2020 Intel Corporation. All rights reserved.
//
// Author: Shriram Shastry <malladi.sastry@linux.intel.com>
#include "typdef.h"
#include <string.h>

/*
 * function x = data_initialization_fixpt
 * Arguments    : int32_t x[TEST_VECTOR]
 * Return Type  : void
 */
 /* 'data_initialization x = fi([-32:0.1:32], 1, 32, 26, fm); Q6.26*/
//Example
//>> tmp = [-2140772762, -2134061876, -2127350989, -2120640103, -2113929216, -2107218330, -2100507444, -2093796557, -2087085671, -2080374784, -2073663898, -2066953012, -2060242125, -2053531239, -2046820352, -2040109466, -2033398580, -2026687693, -2019976807, -2013265920, -2006555034, -1999844148, -1993133261, -1986422375, -1979711488, -1973000602, -1966289716, -1959578829, -1952867943, -1946157056, -1939446170, -1932735284, -1926024397];
//>> tmp / 2 ^ 26
//
//- 31.9000 - 31.8000 - 31.7000 - 31.6000 - 31.5000 - 31.4000 - 31.3000 - 31.2000 - 31.1000 - 31.0000 - 30.9000 - 30.8000 - 30.7000 - 30.6000 - 30.5000 - 30.4000 - 30.3000
//
//- 30.2000 - 30.1000 - 30.0000 - 29.9000 - 29.8000 - 29.7000 - 29.6000 - 29.5000 - 29.4000 - 29.3000 - 29.2000 - 29.1000 - 29.0000 - 28.9000 - 28.8000 - 28.7000

int32_t log10_linear_log[TEST_VECTOR] = {0};
void data_initialization_fixpt(int32_t x[TEST_VECTOR])
{
    static const int32_t iv[TEST_VECTOR] = { -2140772762, -2134061876, -2127350989,
    -2120640103, -2113929216, -2107218330, -2100507444, -2093796557, -2087085671,
    -2080374784, -2073663898, -2066953012, -2060242125, -2053531239, -2046820352,
    -2040109466, -2033398580, -2026687693, -2019976807, -2013265920, -2006555034,
    -1999844148, -1993133261, -1986422375, -1979711488, -1973000602, -1966289716,
    -1959578829, -1952867943, -1946157056, -1939446170, -1932735284, -1926024397,
    -1919313511, -1912602624, -1905891738, -1899180852, -1892469965, -1885759079,
    -1879048192, -1872337306, -1865626420, -1858915533, -1852204647, -1845493760,
    -1838782874, -1832071988, -1825361101, -1818650215, -1811939328, -1805228442,
    -1798517556, -1791806669, -1785095783, -1778384896, -1771674010, -1764963124,
    -1758252237, -1751541351, -1744830464, -1738119578, -1731408692, -1724697805,
    -1717986919, -1711276032, -1704565146, -1697854260, -1691143373, -1684432487,
    -1677721600, -1671010714, -1664299828, -1657588941, -1650878055, -1644167168,
    -1637456282, -1630745396, -1624034509, -1617323623, -1610612736, -1603901850,
    -1597190964, -1590480077, -1583769191, -1577058304, -1570347418, -1563636532,
    -1556925645, -1550214759, -1543503872, -1536792986, -1530082100, -1523371213,
    -1516660327, -1509949440, -1503238554, -1496527668, -1489816781, -1483105895,
    -1476395008, -1469684122, -1462973236, -1456262349, -1449551463, -1442840576,
    -1436129690, -1429418804, -1422707917, -1415997031, -1409286144, -1402575258,
    -1395864372, -1389153485, -1382442599, -1375731712, -1369020826, -1362309940,
    -1355599053, -1348888167, -1342177280, -1335466394, -1328755508, -1322044621,
    -1315333735, -1308622848, -1301911962, -1295201076, -1288490189, -1281779303,
    -1275068416, -1268357530, -1261646644, -1254935757, -1248224871, -1241513984,
    -1234803098, -1228092212, -1221381325, -1214670439, -1207959552, -1201248666,
    -1194537780, -1187826893, -1181116007, -1174405120, -1167694234, -1160983348,
    -1154272461, -1147561575, -1140850688, -1134139802, -1127428916, -1120718029,
    -1114007143, -1107296256, -1100585370, -1093874484, -1087163597, -1080452711,
    -1073741824, -1067030938, -1060320052, -1053609165, -1046898279, -1040187392,
    -1033476506, -1026765620, -1020054733, -1013343847, -1006632960, -999922074,
    -993211188, -986500301, -979789415, -973078528, -966367642, -959656756,
    -952945869, -946234983, -939524096, -932813210, -926102324, -919391437,
    -912680551, -905969664, -899258778, -892547892, -885837005, -879126119,
    -872415232, -865704346, -858993460, -852282573, -845571687, -838860800,
    -832149914, -825439028, -818728141, -812017255, -805306368, -798595482,
    -791884596, -785173709, -778462823, -771751936, -765041050, -758330164,
    -751619277, -744908391, -738197504, -731486618, -724775732, -718064845,
    -711353959, -704643072, -697932186, -691221300, -684510413, -677799527,
    -671088640, -664377754, -657666868, -650955981, -644245095, -637534208,
    -630823322, -624112436, -617401549, -610690663, -603979776, -597268890,
    -590558004, -583847117, -577136231, -570425344, -563714458, -557003572,
    -550292685, -543581799, -536870912, -530160026, -523449140, -516738253,
    -510027367, -503316480, -496605594, -489894708, -483183821, -476472935,
    -469762048, -463051162, -456340276, -449629389, -442918503, -436207616,
    -429496730, -422785844, -416074957, -409364071, -402653184, -395942298,
    -389231412, -382520525, -375809639, -369098752, -362387866, -355676980,
    -348966093, -342255207, -335544320, -328833434, -322122548, -315411661,
    -308700775, -301989888, -295279002, -288568116, -281857229, -275146343,
    -268435456, -261724570, -255013684, -248302797, -241591911, -234881024,
    -228170138, -221459252, -214748365, -208037479, -201326592, -194615706,
    -187904820, -181193933, -174483047, -167772160, -161061274, -154350388,
    -147639501, -140928615, -134217728, -127506842, -120795956, -114085069,
    -107374183, -100663296, -93952410, -87241524, -80530637, -73819751,
    -67108864, -60397978, -53687092, -46976205, -40265319, -33554432, -26843546,
    -20132660, -13421773, -6710887, 0, 6710886, 13421772, 20132659, 26843545,
    33554431, 40265318, 46976204, 53687091, 60397977, 67108863, 73819750,
    80530636, 87241523, 93952409, 100663295, 107374182, 114085068, 120795955,
    127506841, 134217727, 140928614, 147639500, 154350387, 161061273, 167772159,
    174483046, 181193932, 187904819, 194615705, 201326591, 208037478, 214748364,
    221459251, 228170137, 234881023, 241591910, 248302796, 255013683, 261724569,
    268435455, 275146342, 281857228, 288568115, 295279001, 301989887, 308700774,
    315411660, 322122547, 328833433, 335544319, 342255206, 348966092, 355676979,
    362387865, 369098751, 375809638, 382520524, 389231411, 395942297, 402653183,
    409364070, 416074956, 422785843, 429496729, 436207615, 442918502, 449629388,
    456340275, 463051161, 469762047, 476472934, 483183820, 489894707, 496605593,
    503316479, 510027366, 516738252, 523449139, 530160025, 536870911, 543581798,
    550292684, 557003571, 563714457, 570425343, 577136230, 583847116, 590558003,
    597268889, 603979775, 610690662, 617401548, 624112435, 630823321, 637534207,
    644245094, 650955980, 657666867, 664377753, 671088639, 677799526, 684510412,
    691221299, 697932185, 704643071, 711353958, 718064844, 724775731, 731486617,
    738197503, 744908390, 751619276, 758330163, 765041049, 771751935, 778462822,
    785173708, 791884595, 798595481, 805306367, 812017254, 818728140, 825439027,
    832149913, 838860799, 845571686, 852282572, 858993459, 865704345, 872415231,
    879126118, 885837004, 892547891, 899258777, 905969663, 912680550, 919391436,
    926102323, 932813209, 939524095, 946234982, 952945868, 959656755, 966367641,
    973078527, 979789414, 986500300, 993211187, 999922073, 1006632959,
    1013343846, 1020054732, 1026765619, 1033476505, 1040187391, 1046898278,
    1053609164, 1060320051, 1067030937, 1073741823, 1080452710, 1087163596,
    1093874483, 1100585369, 1107296256, 1114007142, 1120718028, 1127428915,
    1134139801, 1140850688, 1147561574, 1154272460, 1160983347, 1167694233,
    1174405120, 1181116006, 1187826892, 1194537779, 1201248665, 1207959552,
    1214670438, 1221381324, 1228092211, 1234803097, 1241513984, 1248224870,
    1254935756, 1261646643, 1268357529, 1275068416, 1281779302, 1288490188,
    1295201075, 1301911961, 1308622848, 1315333734, 1322044620, 1328755507,
    1335466393, 1342177280, 1348888166, 1355599052, 1362309939, 1369020825,
    1375731712, 1382442598, 1389153484, 1395864371, 1402575257, 1409286144,
    1415997030, 1422707916, 1429418803, 1436129689, 1442840576, 1449551462,
    1456262348, 1462973235, 1469684121, 1476395008, 1483105894, 1489816780,
    1496527667, 1503238553, 1509949440, 1516660326, 1523371212, 1530082099,
    1536792985, 1543503872, 1550214758, 1556925644, 1563636531, 1570347417,
    1577058304, 1583769190, 1590480076, 1597190963, 1603901849, 1610612736,
    1617323622, 1624034508, 1630745395, 1637456281, 1644167168, 1650878054,
    1657588940, 1664299827, 1671010713, 1677721600, 1684432486, 1691143372,
    1697854259, 1704565145, 1711276032, 1717986918, 1724697804, 1731408691,
    1738119577, 1744830464, 1751541350, 1758252236, 1764963123, 1771674009,
    1778384896, 1785095782, 1791806668, 1798517555, 1805228441, 1811939328,
    1818650214, 1825361100, 1832071987, 1838782873, 1845493760, 1852204646,
    1858915532, 1865626419, 1872337305, 1879048192, 1885759078, 1892469964,
    1899180851, 1905891737, 1912602624, 1919313510, 1926024396, 1932735283,
    1939446169, 1946157056, 1952867942, 1959578828, 1966289715, 1973000601,
    1979711488, 1986422374, 1993133260, 1999844147, 2006555033, 2013265920,
    2019976806, 2026687692, 2033398579, 2040109465, 2046820352, 2053531238,
    2060242124, 2066953011, 2073663897, 2080374784, 2087085670, 2093796556,
    2100507443, 2107218329, 2113929216, 2120640102, 2127350988, 2134061875,
    2140772761 };

    memcpy(&x[0], &iv[0], 639U * (sizeof(int32_t)));

}
