#' Brazilian Farmers
#'
#' From Valente (1995) \dQuote{In the mid-1960s, Rogers and others conducted
#' an ambitious \sQuote{three country study} to determine influences on adoption
#' of farm practices in Nigeria, India and Brazil. [...] Only in Brazil, and
#' only for hybrid corn, did adoption of the innovation reach more than a small
#' proportion of the farmers.}
#'
#' The dataset has 692 respondents (farmers) from 11 communities. Collected
#' during 1966, it spans 20 years of farming pracitices.
#'
#' @format A data frame with 692 rows and 148 columns:
#' \describe{
#'    \item{village}{village number}
#'    \item{idold}{respondent id}
#'    \item{age}{respondent's age}
#'    \item{liveout}{Lived outside of community}
#'    \item{visits}{# of visits to large city}
#'    \item{contact}{# of contacts with relatives}
#'    \item{coop}{membership in coop}
#'    \item{orgs}{membership in organizations}
#'    \item{patry}{Patriarchalism score}
#'    \item{liter}{Literate}
#'    \item{news1}{# of newspapers or mags pr mon}
#'    \item{subs}{subscribe to news}
#'    \item{radio1}{Own radio}
#'    \item{radio2}{Frequency radio listening}
#'    \item{radio3}{program preference}
#'    \item{tv}{frequency Tv viewing}
#'    \item{movie}{freq movie attendance}
#'    \item{letter}{freq letter writing}
#'    \item{source}{total # of sources used for ag}
#'    \item{practA}{Ever used practice A}
#'    \item{practB}{Ever used practice B}
#'    \item{practC}{Ever used practice C}
#'    \item{practD}{Ever used practice D}
#'    \item{practE}{Ever used practice E}
#'    \item{practF}{Ever used practice F}
#'    \item{practG}{Ever used practice G}
#'    \item{practH}{Ever used practice H}
#'    \item{practI}{Ever used practice I}
#'    \item{practJ}{Ever used practice J}
#'    \item{practK}{Ever used practice K}
#'    \item{practL}{Ever used practice L}
#'    \item{yrA}{A year of adoption}
#'    \item{yrB}{B year of adoption}
#'    \item{yrC}{C year of adoption}
#'    \item{yrD}{D year of adoption}
#'    \item{yrE}{E year of adoption}
#'    \item{yrF}{F year of adoption}
#'    \item{yrG}{G year of adoption}
#'    \item{yrH}{H year of adoption}
#'    \item{yrI}{I year of adoption}
#'    \item{yrJ}{J year of adoption}
#'    \item{yrK}{K year of adoption}
#'    \item{yrL}{L year of adoption}
#'    \item{curA}{A Current use}
#'    \item{curB}{B Current use}
#'    \item{curC}{C Current use}
#'    \item{curD}{D Current use}
#'    \item{curE}{E Current use}
#'    \item{curF}{F Current use}
#'    \item{curG}{G Current use}
#'    \item{curH}{H Current use}
#'    \item{curI}{I Current use}
#'    \item{curJ}{J Current use}
#'    \item{curK}{K Current use}
#'    \item{curL}{L Current use}
#'    \item{srce1}{Source of aware in A}
#'    \item{timeA}{Years ago 1st aware}
#'    \item{src2}{Source of more info on A}
#'    \item{src3}{Most influential source}
#'    \item{use}{use during trial stage}
#'    \item{total}{total # of practices adopted}
#'    \item{futatt}{Future attitude}
#'    \item{achiev}{Achievement Score}
#'    \item{attcred}{Attitude toward credit}
#'    \item{littest}{Score on functional literacy t}
#'    \item{acarcomm}{Communication with ACAR repres}
#'    \item{econk}{Economic knowledge}
#'    \item{caact}{recognize any change agent act}
#'    \item{hfequip}{# of home & farm equips owned}
#'    \item{politk}{political knowledge score}
#'    \item{income}{income}
#'    \item{land1}{total land area in pasture}
#'    \item{land2}{total land area planted}
#'    \item{cows}{# of cows giving milk}
#'    \item{land3}{total land owned}
#'    \item{respf}{respondent named as friend}
#'    \item{respa}{respondent named as ag adv}
#'    \item{resppa}{respondent named for practic A}
#'    \item{resppb}{respondent named for practic B}
#'    \item{resppc}{respondent named for practic C}
#'    \item{poly}{polymorphic OL for 3 practices}
#'    \item{respl}{respondent named for loan}
#'    \item{resppi}{resp named for price info}
#'    \item{repsccp}{resp named for coop comm proj}
#'    \item{counter}{counterfactuality score}
#'    \item{opinion}{opinionness score}
#'    \item{school}{years of schooling by resp}
#'    \item{pk1}{political know 1}
#'    \item{pk2}{political know 2}
#'    \item{pk3}{political know 3}
#'    \item{pk4}{political know 4}
#'    \item{pk5}{political know 5}
#'    \item{innovtim}{innovativeness time}
#'    \item{adoptpct}{adoption percent}
#'    \item{discon}{# of practices discontinued}
#'    \item{mmcred}{Mass media credibility}
#'    \item{trust}{Trust}
#'    \item{stusincn}{Status inconsistency}
#'    \item{nach}{N achievement motivation}
#'    \item{attcred2}{Attitude toward credit}
#'    \item{risk}{Risk taking}
#'    \item{socpart}{Social participate}
#'    \item{patriarc}{patriarchy}
#'    \item{crdit2}{attit to credit for product}
#'    \item{visicit}{visitin cities}
#'    \item{nondep}{non-dependence on farming}
#'    \item{oltotal}{OL total 7 items t-score}
#'    \item{innov}{overall innovativeness score}
#'    \item{icosmo}{cosmo index}
#'    \item{immexp}{mass media exposure index}
#'    \item{iempath}{empathy index}
#'    \item{iach5}{achievement motivation index 5}
#'    \item{iach7}{achievement motivation index 7}
#'    \item{ipk}{political knowledge index}
#'    \item{immc}{mass media credibililty index}
#'    \item{iol}{OL index}
#'    \item{yr}{Actual Year of Adoption}
#'    \item{fs}{ --- MISSING INFO --- }
#'    \item{ado}{Time of Adoption}
#'    \item{tri}{Triangular values used as appro}
#'    \item{hlperc}{high low percent of diffusion}
#'    \item{hlperc1}{ --- MISSING INFO --- }
#'    \item{new}{new or old villages}
#'    \item{card1}{card number}
#'    \item{sour1}{Source: radio}
#'    \item{sour2}{Source: TV}
#'    \item{sour3}{Source: Newpaper}
#'    \item{sour4}{Source: Magazine}
#'    \item{sour5}{Source: ACAR Bulletin}
#'    \item{sour6}{Source: Agronomist}
#'    \item{sour7}{Source: Neighbor}
#'    \item{sourc6}{ --- MISSING INFO --- }
#'    \item{adopt}{ --- MISSING INFO --- }
#'    \item{net31}{nomination friend 1}
#'    \item{net32}{nomination friend 2}
#'    \item{net33}{nomination friend 3}
#'    \item{net21}{nomination influential 1}
#'    \item{net22}{nomination influential 2}
#'    \item{net23}{nomination influential 3}
#'    \item{net11}{nomination practice A}
#'    \item{net12}{nomination practice B}
#'    \item{net13}{nomination practice C}
#'    \item{net41}{nomination coop comm proj}
#'    \item{id}{ --- MISSING INFO --- }
#'    \item{commun}{Number of community}
#'    \item{toa}{Time of Adoption}
#'    \item{test}{ --- MISSING INFO --- }
#'    \item{study}{Number of study in Valente (1995)}
#'  }
#' @source
#' The Brazilian Farmers data were collected as part of a USAID-funded study of farming
#' practicing in the three countries, India, Nigeria, and Brazil.
#' There was only one wave of data that contained survey questions regarding
#' social networks, and only in Brazil did diffusion of the studied farming
#' innovations reach an appreciable saturation level- that was for hybrid seed
#' corn. The data were stored along with hundreds of other datasets by the
#' University of Wisconsin library and I, Tom Valente, paid a fee to have the
#' disks mailed to me in the early 1990s.
#'
#' @references
#'
#' Rogers, E. M., Ascroft, J. R., & Röling, N. (1970). Diffusion of Innovation
#' in Brazil, Nigeria, and India. Unpublished Report. Michigan State University,
#' East Lansing.
#'
#' Valente, T. W. (1995). Network models of the diffusion of innovations (2nd ed.).
#' Cresskill N.J.: Hampton Press.
#' @family diffusion datasets
"brfarmers"

#' \code{diffnet} version of the Brazilian Farmers data
#'
#' A directed dynamic graph with 692 vertices and 21 time periods. The attributes
#' in the graph are static and described in \code{\link{brfarmers}}.
#'
#' @format A \code{\link{diffnet}} class object.
#' @family diffusion datasets
"brfarmersDiffNet"

#' Korean Family Planning
#'
#' From Valente (1995) \dQuote{Scholars at Seoul National University's School
#' of Public Health (Park, Chung, Han & Lee, 1974) collected data on the
#' adoption of family planning methods among all married women of child-bearing
#' age 25 in Korea villages in 1973 (N = 1,047).}
#'
#' The dataset has 1,047 respondents (women) from 25 communities. Collected
#' during 1973 it spans 11 years of data.
#'
#' @format A data frame with 1,047 rows and 432 columns:
#' \describe{
#'    \item{village}{Village of residence}
#'    \item{id}{Respondent ID number}
#'    \item{recno1}{Card number NA}
#'    \item{studno1}{Study number NA}
#'    \item{area1}{Village of residence}
#'    \item{id1}{Respondent ID number}
#'    \item{nmage1}{Number males age 0}
#'    \item{nmage2}{Number males age 0-4}
#'    \item{nmage3}{Number males age 5-9}
#'    \item{nmage4}{Number males age 10-14}
#'    \item{nmage5}{Number males age 15-19}
#'    \item{nmage6}{Number males age 20-24}
#'    \item{nmage7}{Number males age 25-29}
#'    \item{nmage8}{Number males age 30-34}
#'    \item{nmage9}{Number males age 35-39}
#'    \item{nmage10}{Number males age 40-44}
#'    \item{nmage11}{Number males age 45-49}
#'    \item{nmage12}{Number males age 50-54}
#'    \item{nmage13}{Number males age 55-59}
#'    \item{nmage14}{Number males age 60-64}
#'    \item{nmage15}{Number males age 65-69}
#'    \item{nmage16}{Number males age 70-74}
#'    \item{nmage17}{Number males age 75-79}
#'    \item{nmage18}{Number males age 80+}
#'    \item{nfage1}{Number females age 0}
#'    \item{nfage2}{Number females age 0-4}
#'    \item{nfage3}{Number females age 5-9}
#'    \item{nfage4}{Number females age 10-14}
#'    \item{nfage5}{Number females age 15-19}
#'    \item{nfage6}{Number females age 20-24}
#'    \item{nfage7}{Number females age 25-29}
#'    \item{nfage8}{Number females age 30-34}
#'    \item{nfage9}{Number females age 35-39}
#'    \item{nfage10}{Number females age 40-44}
#'    \item{nfage11}{Number females age 45-49}
#'    \item{nfage12}{Number females age 50-54}
#'    \item{nfage13}{Number females age 55-59}
#'    \item{nfage14}{Number females age 60-64}
#'    \item{nfage15}{Number females age 65-69}
#'    \item{nfage16}{Number females age 70-74}
#'    \item{nfage17}{Number females age 75-79}
#'    \item{nfage18}{Number females age 80+}
#'    \item{pregs}{total pregnancies}
#'    \item{pregs1}{number normal deliveries}
#'    \item{pregs2}{number of induced abortions}
#'    \item{pregs3}{number of spontaneous abortions}
#'    \item{pregs4}{number of still births}
#'    \item{pregs5}{number of deaths after live birth}
#'    \item{pregs6}{currently pregnant}
#'    \item{sons}{number of sons}
#'    \item{daughts}{number of daughters}
#'    \item{planning}{Ever heard of FP or birth control}
#'    \item{loop1}{Awareness of Loop}
#'    \item{loop2}{Detailed knowledge of Loop}
#'    \item{loop3}{Attitudes toward Loop}
#'    \item{loop4}{Knowledge of Loop used by neighbors}
#'    \item{loop5}{Knowledge of place of service for Loop}
#'    \item{pill1}{Awareness of Pill}
#'    \item{pill2}{Detailed knowledge of Pill}
#'    \item{pill3}{Attitudes toward Pill}
#'    \item{pill4}{Knowledge of Pill used by neighbors}
#'    \item{pill5}{Knowledge of place of service for Pill}
#'    \item{vase1}{Awareness of Vasectomy}
#'    \item{vase2}{Detailed knowledge of Vasectomy}
#'    \item{vase3}{Attitudes toward Vasectomy}
#'    \item{vase4}{Knowledge of Vasectomy used by neighbors}
#'    \item{vase5}{Knowledge of place of service for Vasectomy}
#'    \item{cond1}{Awareness of Condoms}
#'    \item{cond2}{Detailed knowledge Condoms}
#'    \item{cond3}{Attitudes toward Condoms}
#'    \item{cond4}{Knowledge of Condoms used by neighbors}
#'    \item{cond5}{Knowledge of place of service for Condoms}
#'    \item{rhyt1}{Awareness of Rhythm}
#'    \item{rhyt2}{Detailed knowledge Rhythm}
#'    \item{rhyt3}{Attitudes toward Rhythm}
#'    \item{rhyt4}{Knowledge of Rhythm used by neighbors}
#'    \item{bbt1}{Awareness of Basic Body Temperature}
#'    \item{bbt2}{Detailed knowledge Basic Body Temperature}
#'    \item{bbt3}{Attitudes toward BBT}
#'    \item{recno2}{Record Number NA}
#'    \item{studno2}{Study Number NA}
#'    \item{area2}{village number}
#'    \item{id2}{id number}
#'    \item{bbt4}{Knowledge of BBT used by neighbors}
#'    \item{diap1}{Awareness of Diaphragm}
#'    \item{diap2}{Detailed knowledge Diaphragm}
#'    \item{diap3}{Attitudes toward Diaphragm}
#'    \item{diap4}{Knowledge of Diaphragm used by neighbors}
#'    \item{with1}{Awareness of Withdrawal}
#'    \item{with2}{Detailed knowledge Withdrawal}
#'    \item{with3}{Attitudes toward Withdrawal}
#'    \item{with4}{Knowledge of Withdrawal used by neighbors}
#'    \item{tuba1}{Awareness of Tubal Ligation}
#'    \item{tuba2}{Detailed knowledge TL}
#'    \item{tuba3}{Attitudes toward TL}
#'    \item{tuba4}{Knowledge of TL used by neighbors}
#'    \item{fp1}{Experience with an FP practice}
#'    \item{fp2}{Reasons for not practicing}
#'    \item{fp3}{What would you do if problem was solved}
#'    \item{fp4}{Any other reason for not practicing}
#'    \item{fp5}{Reasons for practicing}
#'    \item{fp6}{time between decision and adoption}
#'    \item{fp7}{reasons for time lag}
#'    \item{fp8}{Ever discontinued practicing}
#'    \item{fp9}{Reasons for discontinuing}
#'    \item{fp10}{Attitude toward FP}
#'    \item{child1}{Ideal number of sons}
#'    \item{child2}{Ideal number of daughters}
#'    \item{child3}{Ideal number of children regardless of sex}
#'    \item{child4}{what do if kept having girls}
#'    \item{comop1}{Spousal communication on # of children}
#'    \item{comop2}{Spousal communication on FP}
#'    \item{comop3}{Consensus on opinion between couple}
#'    \item{comop4}{What was the difference}
#'    \item{comop5}{Opinion on who should practice}
#'    \item{comop6}{Different opinions on who should practice}
#'    \item{comop7}{Who should make final decision}
#'    \item{comop8}{Residence in old age}
#'    \item{net11}{Neighbors talk to about FP- 1}
#'    \item{net12}{Neighbors talk to about FP- 2}
#'    \item{net13}{Neighbors talk to about FP- 3}
#'    \item{net14}{Neighbors talk to about FP- 4}
#'    \item{net15}{Neighbors talk to about FP- 5}
#'    \item{famawe1}{Family members of FP Practice}
#'    \item{famawe2}{Parents awareness of FP Practice}
#'    \item{famawe3}{How did parents-in-law become aware}
#'    \item{famawe4}{How did parents become aware}
#'    \item{famawe5}{How did husband become aware}
#'    \item{advic1}{Advice given to neighbors where to go}
#'    \item{advic2}{Advice given on method}
#'    \item{advic3}{Ever met persons who give advice on FP}
#'    \item{advic4}{Credibility of person advising on FP}
#'    \item{advic5}{Counter advice given to others}
#'    \item{rumor1}{Rumors on Loop}
#'    \item{rumor2}{Rumors on Pill}
#'    \item{rumor3}{Rumors on Vasectomy}
#'    \item{rumor4}{Rumors on Condom}
#'    \item{rumor5}{Rumors on Tuballigation}
#'    \item{media1}{Possession of Radio}
#'    \item{media2}{Possession of TV}
#'    \item{media3}{Subscription to Newspaper}
#'    \item{media4}{Subscription to Happy Home}
#'    \item{media5}{Subscription to other magazine}
#'    \item{media6}{Radio exposure to FP}
#'    \item{media7}{TV exposure to FP}
#'    \item{media8}{Daily paper exposure to FP}
#'    \item{media9}{Happy Home exposure to FP}
#'    \item{media10}{Magazine exposure to FP}
#'    \item{media11}{Movie or slide exposure to FP}
#'    \item{media12}{Poster exposure to FP}
#'    \item{media13}{Pamphlet exposure to FP}
#'    \item{media14}{FP Meeting exposure to FP}
#'    \item{recno3}{Record number NA}
#'    \item{studno3}{Study number NA}
#'    \item{area3}{village}
#'    \item{id3}{id}
#'    \item{media15}{Public lecture exposure to FP}
#'    \item{media16}{Mobile van exposure to FP}
#'    \item{media17}{Neighbors exposure to FP}
#'    \item{media18}{Workers home visiting exposure to FP}
#'    \item{media19}{Husband exposure to FP}
#'    \item{club1}{Awareness of clubs in community}
#'    \item{club2}{Membership in club}
#'    \item{club3}{Reasons for not becoming a member}
#'    \item{club4}{Feeling of necessity of club}
#'    \item{club5}{Visit of mobile van to area}
#'    \item{club6}{Service received from van}
#'    \item{club7}{Decision-making on FP on # children}
#'    \item{club8}{Decision-making on important goods}
#'    \item{club9}{Decision-making on childrens discipline}
#'    \item{club10}{Decision making on purchase wife clothes}
#'    \item{net21}{Closest neighbor most frequently met}
#'    \item{n1adv}{Advice received from neighbor 1}
#'    \item{n1prac}{practice of FP by neighbor 1}
#'    \item{net22}{Closest neighbor person 2}
#'    \item{n2adv}{Advice received from neighbor 2}
#'    \item{n2prac}{Practice of FP by neighbor 2}
#'    \item{net23}{Closest neighbor person 3}
#'    \item{n3adv}{Advice received from neighbor 3}
#'    \item{n3prac}{Practice of FP by neighbor 3}
#'    \item{net24}{Closest neighbor 4}
#'    \item{n4adv}{Advice received from neighbor 4}
#'    \item{n4prac}{Practice of FP by neighbor 4}
#'    \item{net25}{Closest neighbor 5}
#'    \item{n5adv}{Advice received from neighbor 5}
#'    \item{n5prac}{Practice of FP by neighbor 5}
#'    \item{stand}{Standard living of above neighbors}
#'    \item{educ}{Education level of named neighbors}
#'    \item{net31}{Advice on FP sought from 1}
#'    \item{net32}{Advice on FP sought from 2}
#'    \item{net33}{Advice on FP sought from 3}
#'    \item{net34}{Advice on FP sought from 4}
#'    \item{net35}{Advice on FP sought from 5}
#'    \item{net41}{Information provided on FP by 1}
#'    \item{net42}{Information provided on FP by 1}
#'    \item{net43}{Information provided on FP by 1}
#'    \item{net44}{Information provided on FP by 1}
#'    \item{net45}{Information provided on FP by 1}
#'    \item{net51}{Seek advice on induced abortion 1}
#'    \item{net52}{Seek advice on induced abortion 2}
#'    \item{net53}{Seek advice on induced abortion 3}
#'    \item{net54}{Seek advice on induced abortion 4}
#'    \item{net55}{Seek advice on induced abortion 5}
#'    \item{age}{Age of respondent}
#'    \item{agemar}{Age at first marriage}
#'    \item{recno4}{Rec no NA}
#'    \item{studno4}{Study no NA}
#'    \item{area4}{village}
#'    \item{id4}{id}
#'    \item{net61}{Advice on health sought from 1}
#'    \item{net62}{Advice on health sought from 2}
#'    \item{net63}{Advice on health sought from 3}
#'    \item{net64}{Advice on health sought from 4}
#'    \item{net65}{Advice on health sought from 5}
#'    \item{net71}{Advice on purchase of goods 1}
#'    \item{net72}{Advice on purchase of goods 2}
#'    \item{net73}{Advice on purchase of goods 3}
#'    \item{net74}{Advice on purchase of goods 4}
#'    \item{net75}{Advice on purchase of goods 5}
#'    \item{net81}{Advice on childrens education 1}
#'    \item{net82}{Advice on childrens education 2}
#'    \item{net83}{Advice on childrens education 3}
#'    \item{net84}{Advice on childrens education 4}
#'    \item{net85}{Advice on childrens education 5}
#'    \item{rfampl1}{Advice on FP sought by 1}
#'    \item{rfampl2}{Advice on FP sought by 2}
#'    \item{rfampl3}{Advice on FP sought by 3}
#'    \item{rfampl4}{Advice on FP sought by 4}
#'    \item{rfampl5}{Advice on FP sought by 5}
#'    \item{rfampll}{Leadership score - indegree FP}
#'    \item{rabort1}{Advice on abortion sought by 1}
#'    \item{rabort2}{Advice on abortion sought by 2}
#'    \item{rabort3}{Advice on abortion sought by 3}
#'    \item{rabort4}{Advice on abortion sought by 4}
#'    \item{rabort5}{Advice on abortion sought by 5}
#'    \item{rabortl}{Leadership score - indegree abortion}
#'    \item{rhealth1}{Advice on health sought by 1}
#'    \item{rhealth2}{Advice on health sought by}
#'    \item{rhealth3}{Advice on health sought by}
#'    \item{rhealth4}{Advice on health sought by}
#'    \item{rhealth5}{Advice on health sought by}
#'    \item{rhealthl}{Leadership score - indegree health}
#'    \item{recno5}{rec no NA}
#'    \item{studno5}{study no NA}
#'    \item{area5}{village}
#'    \item{id5}{id}
#'    \item{rgoods1}{Advice on purchases sought by 1}
#'    \item{rgoods2}{Advice on purchases sought by 2}
#'    \item{rgoods3}{Advice on purchases sought by 3}
#'    \item{rgoods4}{Advice on purchases sought by 4}
#'    \item{rgoods5}{Advice on purchases sought by 5}
#'    \item{rgoodsl}{Leadership score - indegree purchases}
#'    \item{reduc1}{Advice on education sought by 1}
#'    \item{reduc2}{Advice on education sought by 2}
#'    \item{reduc3}{Advice on education sought by 3}
#'    \item{reduc4}{Advice on education sought by 4}
#'    \item{reduc5}{Advice on education sought by 5}
#'    \item{reducl}{Leadership score - indegree education}
#'    \item{hub1}{Husbands friend 1}
#'    \item{hub2}{Husbands friend 2}
#'    \item{hub3}{Husbands friend 3}
#'    \item{hub4}{Husbands friend 4}
#'    \item{hub5}{Husbands friend 5}
#'    \item{hubed}{Husbands education}
#'    \item{wifeed}{Wifes education}
#'    \item{wiferel}{Wifes religion}
#'    \item{hubocc}{Husbands occupation}
#'    \item{wifeocc}{Wifes occupation}
#'    \item{know1}{Can you insert a loop yourself}
#'    \item{know2}{Can you remove it alone}
#'    \item{know3}{Can a man use a loop}
#'    \item{know4}{How long can a loop be used}
#'    \item{know5}{Which doctor}
#'    \item{know6}{Doctor or nurse}
#'    \item{know7}{Oral pill method}
#'    \item{know8}{Can men take pills}
#'    \item{know9}{Long term use}
#'    \item{know10}{Time required for vasectomy}
#'    \item{know11}{Does vasectomy = castration}
#'    \item{know12}{Can any doctor do vasectomies}
#'    \item{pref1}{Who prefer use: Husband or wife}
#'    \item{pref2}{Reasons for preferring FP practice by wife}
#'    \item{pref3}{Reasons for preferring FP practice by husband}
#'    \item{ageend}{Ideal age to end childbearing}
#'    \item{cfp}{Current status of FP}
#'    \item{cfatt1}{Husbands attitude}
#'    \item{cfatt2}{In-laws attitude}
#'    \item{cfatt3}{Own parents attitude}
#'    \item{cbyr}{Start of period from year}
#'    \item{cbmnth}{Start of period from month}
#'    \item{ceyr}{End of period year}
#'    \item{cemnth}{End of period month}
#'    \item{clngth}{Length of period}
#'    \item{cawe1}{FP contact}
#'    \item{cawe2}{Awareness of contraceptive method at the time}
#'    \item{cawe3}{Awareness of service site}
#'    \item{cawe4}{Credibiilty}
#'    \item{recno6}{rec no NA}
#'    \item{studno6}{study no NA}
#'    \item{area6}{village}
#'    \item{id6}{id}
#'    \item{fpt1}{FP Status time 1}
#'    \item{fatt1t1}{Husbands attitude T1}
#'    \item{fatt2t1}{In-laws attitude T1}
#'    \item{fatt3t1}{Own parents attitude T1}
#'    \item{byrt1}{Start of Time 1 from year}
#'    \item{lngtht1}{Length of Time 1}
#'    \item{awe1t1}{FP Contact Time 1}
#'    \item{awe2t1}{Methods known at Time 1}
#'    \item{awe3t1}{Knowledge of service sites Time 1}
#'    \item{awe4t1}{Credibility of service site Time 1}
#'    \item{fpt2}{FP Status time 2}
#'    \item{fatt1t2}{Husbands attitude T2}
#'    \item{fatt2t2}{In-laws attitude T2}
#'    \item{fatt3t2}{Own parents attitude T2}
#'    \item{byrt2}{Start of Time 2 from year}
#'    \item{lngtht2}{Length of Time 2}
#'    \item{awe1t2}{FP Contact Time 2}
#'    \item{awe2t2}{Methods known at Time 2}
#'    \item{awe3t2}{Knowledge of service sites Time 2}
#'    \item{awe4t2}{Credibility of service site Time 2}
#'    \item{fpt3}{FP Status time 3}
#'    \item{fatt1t3}{Husbands attitude T3}
#'    \item{fatt2t3}{In-laws attitude T3}
#'    \item{fatt3t3}{Own parents attitude T3}
#'    \item{byrt3}{Start of Time 3 from year}
#'    \item{lngtht3}{Length of Time 3}
#'    \item{awe1t3}{FP Contact Time 3}
#'    \item{awe2t3}{Methods known at Time 3}
#'    \item{awe3t3}{Knowledge of service sites Time 3}
#'    \item{awe4t3}{Credibility of service site Time 3}
#'    \item{fpt4}{FP Status time 4}
#'    \item{fatt1t4}{Husbands attitude T4}
#'    \item{fatt2t4}{In-laws attitude T4}
#'    \item{fatt3t4}{Own parents attitude T4}
#'    \item{byrt4}{Start of Time 4 from year}
#'    \item{lngtht4}{Length of Time 4}
#'    \item{awe1t4}{FP Contact Time 4}
#'    \item{awe2t4}{Methods known at Time 4}
#'    \item{awe3t4}{Knowledge of service sites Time 4}
#'    \item{awe4t4}{Credibility of service site Time 4}
#'    \item{fpt5}{FP Status time 5}
#'    \item{fatt1t5}{Husbands attitude T5}
#'    \item{fatt2t5}{In-laws attitude T5}
#'    \item{fatt3t5}{Own parents attitude T5}
#'    \item{byrt5}{Start of Time 5 from year}
#'    \item{lngtht5}{Length of Time 5}
#'    \item{awe1t5}{FP Contact Time 5}
#'    \item{awe2t5}{Methods known at Time 5}
#'    \item{awe3t5}{Knowledge of service sites Time 5}
#'    \item{awe4t5}{Credibility of service site Time 5}
#'    \item{fpt6}{FP Status time 6}
#'    \item{fatt1t6}{Husbands attitude T6}
#'    \item{fatt2t6}{In-laws attitude T6}
#'    \item{fatt3t6}{Own parents attitude T6}
#'    \item{byrt6}{Start of Time 6 from year}
#'    \item{lngtht6}{Length of Time 6}
#'    \item{awe1t6}{FP Contact Time 6}
#'    \item{awe2t6}{Methods known at Time 6}
#'    \item{awe3t6}{Knowledge of service sites Time 6}
#'    \item{awe4t6}{Credibility of service site Time 6}
#'    \item{recno7}{rec no NA}
#'    \item{studno7}{study no NA}
#'    \item{area7}{village}
#'    \item{id7}{id}
#'    \item{fpt7}{FP Status time 7}
#'    \item{fatt1t7}{Husbands attitude T7}
#'    \item{fatt2t7}{In-laws attitude T7}
#'    \item{fatt3t7}{Own parents attitude T7}
#'    \item{byrt7}{Start of Time 7 from year}
#'    \item{lngtht7}{Length of Time 7}
#'    \item{awe1t7}{FP Contact Time 7}
#'    \item{awe2t7}{Methods known at Time 7}
#'    \item{awe3t7}{Knowledge of service sites Time 7}
#'    \item{awe4t7}{Credibility of service site Time 7}
#'    \item{fpt8}{FP Status time 8}
#'    \item{fatt1t8}{Husbands attitude T8}
#'    \item{fatt2t8}{In-laws attitude T8}
#'    \item{fatt3t8}{Own parents attitude T8}
#'    \item{byrt8}{Start of Time 8 from year}
#'    \item{lngtht8}{Length of Time 8}
#'    \item{awe1t8}{FP Contact Time 8}
#'    \item{awe2t8}{Methods known at Time 8}
#'    \item{awe3t8}{Knowledge of service sites Time 8}
#'    \item{awe4t8}{Credibility of service site Time 8}
#'    \item{fpt9}{FP Status time 9}
#'    \item{fatt1t9}{Husbands attitude T9}
#'    \item{fatt2t9}{In-laws attitude T9}
#'    \item{fatt3t9}{Own parents attitude T9}
#'    \item{byrt9}{Start of Time 9 from year}
#'    \item{lngtht9}{Length of Time 9}
#'    \item{awe1t9}{FP Contact Time 9}
#'    \item{awe2t9}{Methods known at Time 9}
#'    \item{awe3t9}{Knowledge of service sites Time 9}
#'    \item{awe4t9}{Credibility of service site Time 9}
#'    \item{fpt10}{FP Status time 10}
#'    \item{fatt1t10}{Husbands attitude T10}
#'    \item{fatt2t10}{In-laws attitude T10}
#'    \item{fatt3t10}{Own parents attitude T10}
#'    \item{byrt10}{Start of Time 10 from year}
#'    \item{lngtht10}{Length of Time 10}
#'    \item{awe1t10}{FP Contact Time 10}
#'    \item{awe2t10}{Methods known at Time 10}
#'    \item{awe3t10}{Knowledge of service sites Time 10}
#'    \item{awe4t10}{Credibility of service site Time 10}
#'    \item{fpt11}{FP Status time 11}
#'    \item{fatt1t11}{Husbands attitude T11}
#'    \item{fatt2t11}{In-laws attitude T11}
#'    \item{fatt3t11}{Own parents attitude T11}
#'    \item{byrt11}{Start of Time 11 from year}
#'    \item{lngtht11}{Length of Time 11}
#'    \item{awe1t11}{FP Contact Time 11}
#'    \item{awe2t11}{Methods known at Time 11}
#'    \item{awe3t11}{Knowledge of service sites Time 11}
#'    \item{awe4t11}{Credibility of service site Time 11}
#'    \item{fpt12}{FP Status time 12}
#'    \item{fatt1t12}{Husbands attitude T12}
#'    \item{fatt2t12}{In-laws attitude T12}
#'    \item{fatt3t12}{Own parents attitude T12}
#'    \item{byrt12}{Start of Time 12 from year}
#'    \item{lngtht12}{Length of Time 12}
#'    \item{awe1t12}{FP Contact Time 12}
#'    \item{awe2t12}{Methods known at Time 12}
#'    \item{awe3t12}{Knowledge of service sites Time 12}
#'    \item{awe4t12}{Credibility of service site Time 12}
#'    \item{ado}{adopt times years converted to 1=63}
#'    \item{ado1}{}
#'    \item{ado2}{}
#'    \item{ado3}{}
#'    \item{commun}{Village number}
#'    \item{toa}{Time of Adoption}
#'    \item{study}{Study (for when multiple diff studies used)}
#' }
#'
#' @source
#' The Korean Family Planning data were
#' stored on a Vax tape that Rogers had given to Marc Granovetter who then gave
#' it to his colleague Roland Soong (see Granovetter & Soong, 1983).  Granovetter
#' instructed Song to send the tape to me and I had it loaded on the Vax machine
#' at USC in 1990 and was able to download the data to a PC. The first two datasets
#' were acquired for my dissertation (Valente, 1991) and the third added as I
#' completed my book on Network Models of the Diffusion of Innovations (Valente,
#' 1995; also see Valente, 2005).
#' @references
#' Everett M. Rogers, & Kincaid, D. L. (1981). Communication Networks: Toward a
#' New Paradigm for Research. (C. Macmillan, Ed.). New York; London: Free Press.
#'
#' Valente, T. W. (1995). Network models of the diffusion of innovations (2nd ed.).
#' Cresskill N.J.: Hampton Press.
#' @family diffusion datasets
"kfamily"

# "kfamilyDiffNet"

#' \code{diffnet} version of the Korean Family Planning data
#'
#' A directed dynamic graph with 1,047 vertices and 11 time periods. The attributes
#' in the graph are static and described in \code{\link{kfamily}}.
#'
#' @format A \code{\link{diffnet}} class object.
#' @family diffusion datasets
"kfamilyDiffNet"



#' Medical Innovation
#'
#' From Valente (1995) \dQuote{Coleman, Katz and Menzel from Columbia University's Bureau of Applied Research
#' studied the adoption of tetracycline by physiciams in four Illinois communities
#' in 1954.[...] Tetracycline was a powerful and useful antibiotic just introduced in
#' the mid-1950s}
#'
#'
#' The collected dataset has 125 respondents (doctors), and spans 17 months of data
#' collected in 1955. Time of adoption of non-adopters has been set to month
#' 18 (see the manual entry titled \code{\link[netdiffuseR:diffusion-data]{Difussion Network Datasets}}).
#'
#' @format A data frame with 125 rows and 59 columns:
#' \describe{
#'    \item{city}{city id}
#'    \item{id}{sequential respondent id}
#'    \item{detail}{detail man}
#'    \item{meet}{meetings, lectures, hospitals}
#'    \item{coll}{colleagues}
#'    \item{attend}{attend professional meets}
#'    \item{proage}{professional age}
#'    \item{length}{lenght of reside in community}
#'    \item{here}{only practice here}
#'    \item{science}{science versus patients}
#'    \item{position}{position in home base}
#'    \item{journ2}{journal subscriptions}
#'    \item{paadico}{Percent alter adoption date imp}
#'    \item{ado}{adoption month 1 to 18}
#'    \item{thresh}{threshold}
#'    \item{ctl}{corrected tl tl-exp level}
#'    \item{catbak}{category 1-init 2-marg 3-low tl}
#'    \item{sourinfo}{source of information}
#'    \item{origid}{original respondent id}
#'    \item{adopt}{adoption date 1= 11/53}
#'    \item{recon}{reconstructed med innov}
#'    \item{date}{date became aware}
#'    \item{info}{information source}
#'    \item{most}{most important info source}
#'    \item{journ}{journals}
#'    \item{drug}{drug houses}
#'    \item{net1_1}{advisor nomination1}
#'    \item{net1_2}{advisor nomination2}
#'    \item{net1_3}{advisor nomination3}
#'    \item{net2_1}{discuss nomination1}
#'    \item{net2_2}{discuss nomination2}
#'    \item{net2_3}{discuss nomination3}
#'    \item{net3_1}{friends nomination1}
#'    \item{net3_2}{friends nomination2}
#'    \item{net3_3}{friends nomination3}
#'    \item{nojourn}{number of pro journals receive}
#'    \item{free}{free time companions}
#'    \item{social}{med discussions during social}
#'    \item{club}{club membership}
#'    \item{friends}{friends are doctors}
#'    \item{young}{young patients}
#'    \item{nonpoor}{nonpoverty patients}
#'    \item{office}{office visits}
#'    \item{house}{house calls}
#'    \item{tend}{tendency to prescribe drugs}
#'    \item{reltend}{relative tendency to prescribe}
#'    \item{perc}{perceived drug competition}
#'    \item{proximty}{physical proximity to other doc}
#'    \item{home}{home base hospital affiliation}
#'    \item{special}{specialty}
#'    \item{belief}{belief in science}
#'    \item{proage2}{profesional age 2}
#'    \item{presc}{prescription prone}
#'    \item{detail2}{contact with detail man}
#'    \item{dichot}{dichotomous personal preference}
#'    \item{expect}{adoption month expected}
#'    \item{recall}{recalls adopting}
#'    \item{commun}{Number of community}
#'    \item{toa}{Time of Adoption}
#'    \item{study}{Number of study in Valente (1995)}
#' }
#' @source
#' The Medical Innovation data were stored in file cabinets in a basement
#' building at Columbia University. Ron Burt (1987) acquired an NSF grant to
#' develop network diffusion models and retrieve the original surveys and enter
#' them into a database. He distributed copies of the data on diskette and sent
#' one to me, Tom Valente, and I imported onto a PC environment.
#' @references
#' Coleman, J., Katz, E., & Menzel, H. (1966). Medical innovation: A diffusion
#' study (2nd ed.). New York: Bobbs-Merrill
#'
#' Valente, T. W. (1995). Network models of the diffusion of innovations (2nd ed.).
#' Cresskill N.J.: Hampton Press.
#' @family diffusion datasets
"medInnovations"


#' \code{diffnet} version of the Medical Innovation data
#'
#' A directed dynamic graph with 125 vertices and 18 time periods. The attributes
#' in the graph are static and described in \code{\link{medInnovations}}.
#'
#' @format A \code{\link{diffnet}} class object.
#' @family diffusion datasets
"medInnovationsDiffNet"

#' Diffusion Network Datasets
#'
#' @details
#' The three classic network diffusion datasets included in netdiffuseR are the
#' medical innovation data originally collected by Coleman, Katz & Menzel (1966);
#' the Brazilian Farmers collected as part of the three country study implemented
#' by Everett Rogers (Rogers, Ascroft, & Röling, 1970), and Korean Family Planning
#' data collected by researchers at the Seoul National University's School of
#' Public (Rogers & Kincaid, 1981). The table below summarizes the three datasets:
#'
#' \tabular{lccc}{
#'		\tab	\bold{Medical Innovation}	\tab	\bold{Brazilian Farmers}	\tab	\bold{Korean Family Planning}	\cr
#'	\emph{Country}	\tab	USA	\tab	Brazil	\tab	Korean	\cr
#'	\emph{# Respondents}	\tab	125 Doctors	\tab	692 Farmers	\tab	1,047 Women	\cr
#'	\emph{# Communities}	\tab	4	\tab	11	\tab	25	\cr
#'	\emph{Innovation}	\tab	Tetracycline	\tab	Hybrid Corn Seed	\tab	Family Planning	\cr
#'	\emph{Time for Diffusion}	\tab	18 Months	\tab	20 Years	\tab	11 Years	\cr
#'	\emph{Year Data Collected}	\tab	1955-1956	\tab	1966	\tab	1973	\cr
#'	\emph{Ave. Time to 50\%}	\tab	6	\tab	16	\tab	7	\cr
#'	\emph{Highest Saturation}	\tab	0.89	\tab	0.98	\tab	0.83	\cr
#'	\emph{Lowest Saturation}	\tab	0.81	\tab	0.29	\tab	0.44	\cr
#'	\emph{Citation}	\tab	Coleman et al (1966)	\tab	Rogers et al (1970)	\tab	Rogers & Kincaid (1981)	\cr
#'	}
#'
#' All datasets include a column called \emph{study} which is coded as
#' (1) Medical Innovation (2) Brazilian Farmers, (3) Korean Family Planning.
#'
#' @section Right censored data:
#' By convention, non-adopting actors are coded as one plus the last observed time
#' of adoption.  Prior empirical event history approaches have used this approach
#' (Valente, 2005; Marsden and Podolny, 1990) and studies have shown that
#' omitting such observations leads to biased results (van den Bulte & Iyengar,
#' 2011).
#'
#' @references
#' Burt, R. S. (1987). "Social Contagion and Innovation: Cohesion versus
#' Structural Equivalence". American Journal of Sociology, 92(6), 1287–1335.
#' \url{http://doi.org/10.1086/228667}
#'
#' Coleman, J., Katz, E., & Menzel, H. (1966). Medical innovation: A diffusion
#' study (2nd ed.). New York: Bobbs-Merrill
#'
#' Granovetter, M., & Soong, R. (1983). Threshold models of diffusion and
#' collective behavior. The Journal of Mathematical Sociology, 9(October 2013),
#' 165–179. \url{http://doi.org/10.1080/0022250X.1983.9989941}
#'
#' Rogers, E. M., Ascroft, J. R., & Röling, N. (1970). Diffusion of Innovation
#' in Brazil, Nigeria, and India. Unpublished Report. Michigan State University,
#' East Lansing.
#'
#' Everett M. Rogers, & Kincaid, D. L. (1981). Communication Networks: Toward a
#' New Paradigm for Research. (C. Macmillan, Ed.). New York; London: Free Press.
#'
#' Mardsen, P., & Podolny, J. (1990). Dynamic Analysis of Network Diffusion Processes,
#' J. Weesie, H. Flap, eds. Social Networks Through Time, 197–214.
#'
#' Marsden, P. V., & Friedkin, N. E. (1993). Network Studies of Social Influence.
#' Sociological Methods & Research, 22(1), 127–151.
#' \url{http://doi.org/10.1177/0049124193022001006}
#'
#' Van den Bulte, C., & Iyengar, R. (2011). Tricked by Truncation: Spurious
#' Duration Dependence and Social Contagion in Hazard Models. Marketing Science,
#' 30(2), 233–248. \url{http://doi.org/10.1287/mksc.1100.0615}
#'
#' Valente, T. W. (1991). Thresholds and the critical mass: Mathematical models
#' of the diffusion of innovations. University of Southern California.
#'
#' Valente, T. W. (1995). "Network models of the diffusion of innovations" (2nd ed.).
#' Cresskill N.J.: Hampton Press.
#'
#' Valente, T. W. (2005). Network Models and Methods for Studying the Diffusion of Innovations.
#' In Models and Methods in Social Network Analysis, Volume 28 of Structural
#' Analysis in the Social Sciences (pp. 98–116). New York: Cambridge University Press.
#' @family diffusion datasets
#' @name diffusion-data
#' @author Thomas W. Valente
NULL


#' Fake survey data
#'
#' This data frame is used to ilustrate some of the functions of the package,
#' in particular, the \code{\link{survey_to_diffnet}} function. This dataset
#' can be merged with the \code{\link{fakeEdgelist}}.
#'
#' @format A data frame with 9 rows and 9 variables
#' \describe{
#'  \item{id}{Unique id at group level}
#'  \item{toa}{Time of adoption}
#'  \item{group}{Group id}
#'  \item{net1}{Network nomination 1}
#'  \item{net2}{Network nomination 2}
#'  \item{net3}{Network nomination 3}
#'  \item{age}{Age of the respondent}
#'  \item{gender}{Gende of the respondent}
#'  \item{note}{Descroption of the respondent}
#' }
#'
#' @source Generated for the package.
#' @family diffusion datasets
#' @author George G. Vega Yon
"fakesurvey"

#' Fake longitudinal survey data
#'
#' This data frame is used to ilustrate some of the functions of the package,
#' in particular, the \code{\link{survey_to_diffnet}} function. This dataset
#' can be merged with the \code{\link{fakeDynEdgelist}}.
#'
#' @format A data frame with 18 rows and 10 variables
#' \describe{
#'  \item{id}{Unique id at group level}
#'  \item{toa}{Time of adoption}
#'  \item{group}{Group id}
#'  \item{net1}{Network nomination 1}
#'  \item{net2}{Network nomination 2}
#'  \item{net3}{Network nomination 3}
#'  \item{age}{Age of the respondent}
#'  \item{gender}{Gende of the respondent}
#'  \item{note}{Descroption of the respondent}
#'  \item{time}{Timing of the wave}
#' }
#'
#' @source Generated for the package.
#' @family diffusion datasets
#' @author George G. Vega Yon
"fakesurveyDyn"


#' Fake dynamic edgelist
#'
#' A data frame used for examples in reading edgelist format networks. This
#' edgelist can be merged with the dataset \code{\link{fakesurveyDyn}}.
#'
#' @format A data frame with 22 rows and 4 variables
#' \describe{
#'  \item{ego}{Nominating individual}
#'  \item{alter}{Nominated individual}
#'  \item{value}{Strength of the tie}
#'  \item{time}{Integer with the time of the spell}
#' }
#'
#' @source Generated for the package
#' @family diffusion datasets
#' @author George G. Vega Yon
"fakeDynEdgelist"

#' Fake static edgelist
#'
#' A data frame used for examples in reading edgelist format networks. This
#' edgelist can be merged with the dataset \code{\link{fakesurvey}}.
#'
#' @format A data frame with 11 rows and 3 variables
#' \describe{
#'  \item{ego}{Nominating individual}
#'  \item{alter}{Nominated individual}
#'  \item{value}{Strength of the tie}
#' }
#'
#' @source Generated for the package
#' @family diffusion datasets
#' @author George G. Vega Yon
"fakeEdgelist"
