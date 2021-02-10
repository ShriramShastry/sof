// SPDX - License - Identifier: BSD - 3 - Clause
//
//Copyright(c) 2020 Intel Corporation.All rights reserved.
//
//Author : Shriram Shastry <malladi.sastry@linux.intel.com>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
/* Function Definitions */
/*
 * function [x,y]= init_data_fixpt
 * Arguments    : uint32_T x[320]
 *                int32_T y[41]
 * Return Type  : void
 */
void init_data_fixpt(uint32_t x[320], int32_t y[41])
{
    static const int32_t iv[41] = { -1073741824, -1020054733, -966367642,
      -912680550, -858993459, -805306368, -751619277, -697932186, -644245094,
      -590558003, -536870912, -483183821, -429496730, -375809638, -322122547,
      -268435456, -214748365, -161061274, -107374182, -53687091, 0, 53687091,
      107374182, 161061274, 214748365, 268435456, 322122547, 375809638, 429496730,
      483183821, 536870912, 590558003, 644245094, 697932186, 751619277, 805306368,
      858993459, 912680550, 966367642, 1020054733, 1073741824 };

    static const uint32_t uv[320] = { 6710886U, 13421773U, 20132659U, 26843546U,
      33554432U, 40265318U, 46976205U, 53687091U, 60397978U, 67108864U, 73819750U,
      80530637U, 87241523U, 93952410U, 100663296U, 107374182U, 114085069U,
      120795955U, 127506842U, 134217728U, 140928614U, 147639501U, 154350387U,
      161061274U, 167772160U, 174483046U, 181193933U, 187904819U, 194615706U,
      201326592U, 208037478U, 214748365U, 221459251U, 228170138U, 234881024U,
      241591910U, 248302797U, 255013683U, 261724570U, 268435456U, 275146342U,
      281857229U, 288568115U, 295279002U, 301989888U, 308700774U, 315411661U,
      322122547U, 328833434U, 335544320U, 342255206U, 348966093U, 355676979U,
      362387866U, 369098752U, 375809638U, 382520525U, 389231411U, 395942298U,
      402653184U, 409364070U, 416074957U, 422785843U, 429496730U, 436207616U,
      442918502U, 449629389U, 456340275U, 463051162U, 469762048U, 476472934U,
      483183821U, 489894707U, 496605594U, 503316480U, 510027366U, 516738253U,
      523449139U, 530160026U, 536870912U, 543581798U, 550292685U, 557003571U,
      563714458U, 570425344U, 577136230U, 583847117U, 590558003U, 597268890U,
      603979776U, 610690662U, 617401549U, 624112435U, 630823322U, 637534208U,
      644245094U, 650955981U, 657666867U, 664377754U, 671088640U, 677799526U,
      684510413U, 691221299U, 697932186U, 704643072U, 711353958U, 718064845U,
      724775731U, 731486618U, 738197504U, 744908390U, 751619277U, 758330163U,
      765041050U, 771751936U, 778462822U, 785173709U, 791884595U, 798595482U,
      805306368U, 812017254U, 818728141U, 825439027U, 832149914U, 838860800U,
      845571686U, 852282573U, 858993459U, 865704346U, 872415232U, 879126118U,
      885837005U, 892547891U, 899258778U, 905969664U, 912680550U, 919391437U,
      926102323U, 932813210U, 939524096U, 946234982U, 952945869U, 959656755U,
      966367642U, 973078528U, 979789414U, 986500301U, 993211187U, 999922074U,
      1006632960U, 1013343846U, 1020054733U, 1026765619U, 1033476506U, 1040187392U,
      1046898278U, 1053609165U, 1060320051U, 1067030938U, 1073741824U, 1080452710U,
      1087163597U, 1093874483U, 1100585370U, 1107296256U, 1114007142U, 1120718029U,
      1127428915U, 1134139802U, 1140850688U, 1147561574U, 1154272461U, 1160983347U,
      1167694234U, 1174405120U, 1181116006U, 1187826893U, 1194537779U, 1201248666U,
      1207959552U, 1214670438U, 1221381325U, 1228092211U, 1234803098U, 1241513984U,
      1248224870U, 1254935757U, 1261646643U, 1268357530U, 1275068416U, 1281779302U,
      1288490189U, 1295201075U, 1301911962U, 1308622848U, 1315333734U, 1322044621U,
      1328755507U, 1335466394U, 1342177280U, 1348888166U, 1355599053U, 1362309939U,
      1369020826U, 1375731712U, 1382442598U, 1389153485U, 1395864371U, 1402575258U,
      1409286144U, 1415997030U, 1422707917U, 1429418803U, 1436129690U, 1442840576U,
      1449551462U, 1456262349U, 1462973235U, 1469684122U, 1476395008U, 1483105894U,
      1489816781U, 1496527667U, 1503238554U, 1509949440U, 1516660326U, 1523371213U,
      1530082099U, 1536792986U, 1543503872U, 1550214758U, 1556925645U, 1563636531U,
      1570347418U, 1577058304U, 1583769190U, 1590480077U, 1597190963U, 1603901850U,
      1610612736U, 1617323622U, 1624034509U, 1630745395U, 1637456282U, 1644167168U,
      1650878054U, 1657588941U, 1664299827U, 1671010714U, 1677721600U, 1684432486U,
      1691143373U, 1697854259U, 1704565146U, 1711276032U, 1717986918U, 1724697805U,
      1731408691U, 1738119578U, 1744830464U, 1751541350U, 1758252237U, 1764963123U,
      1771674010U, 1778384896U, 1785095782U, 1791806669U, 1798517555U, 1805228442U,
      1811939328U, 1818650214U, 1825361101U, 1832071987U, 1838782874U, 1845493760U,
      1852204646U, 1858915533U, 1865626419U, 1872337306U, 1879048192U, 1885759078U,
      1892469965U, 1899180851U, 1905891738U, 1912602624U, 1919313510U, 1926024397U,
      1932735283U, 1939446170U, 1946157056U, 1952867942U, 1959578829U, 1966289715U,
      1973000602U, 1979711488U, 1986422374U, 1993133261U, 1999844147U, 2006555034U,
      2013265920U, 2019976806U, 2026687693U, 2033398579U, 2040109466U, 2046820352U,
      2053531238U, 2060242125U, 2066953011U, 2073663898U, 2080374784U, 2087085670U,
      2093796557U, 2100507443U, 2107218330U, 2113929216U, 2120640102U, 2127350989U,
      2134061875U, 2140772762U, 2147483648U };

    memcpy(&x[0], &uv[0], 320U * (sizeof(uint32_t)));
    memcpy(&y[0], &iv[0], 41U * (sizeof(int32_t)));
 }
