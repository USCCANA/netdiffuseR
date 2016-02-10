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
#'    \item{commun}{ --- MISSING INFO --- }
#'    \item{toa}{ --- MISSING INFO --- }
#'    \item{test}{ --- MISSING INFO --- }
#'    \item{study}{ --- MISSING INFO --- }
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
#' age 25 in Korea villages in 1973 (N = 1,047)}
#'
#' The dataset has 1,047 respondents (women) from 25 communities. Collected
#' during 1973 it spans 11 years of data.
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
#' @name kp
NULL

# "kfamiliesDiffNet"

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
#'    \item{commun}{}
#'    \item{toa}{}
#'    \item{study}{}
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
#' @author Valente, T.W.
NULL
