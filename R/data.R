#' Brazilian Farmers
#'
#' 692 respondents farmers on 11 communities. Data on the adoption of Hybrid corn
#' with a time span of 20 years. Collected during 1966
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
#' @source SOURCE?
#' @references
#'
#' Rogers, E. M., Ascroft, J. R., & Röling, N. (1970). Diffusion of Innovation
#' in Brazil, Nigeria, and India. Unpublished Report. Michigan State University,
#' East Lansing.
#'
#' Valente, T. W. (1995). Network models of the diffusion of innovations (2nd ed.).
#' Cresskill N.J.: Hampton Press.
"brfarmers"

#' \code{diffnet} version of the Brazilian Farmers data
#'
#' A directed dynamic graph with 692 vertices and 21 time periods. The attributes
#' in the graph are static and described in \code{\link{brfarmers}}.
#'
#' @format A \code{\link{diffnet}} class object.
"brfarmersDiffNet"

#' Korean Family Planning
#'
#' 1047 respondents women
#' 25 communities
#' Family Planning
#' Length of diffusion 11 years
#' Year of collection 1973
#' Years span
#'
#' @source Rogers & Kincaid 1981
#' @references TOM
#' @name kp
NULL

# "kfamiliesDiffNet"

#' Medical Innovation
#'
#' 125 respondents doctors. 4 communities. Innovation: Adoption of Tetracycline.
#' Time span of 18 months collected in 1955.
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
#' @source Coleman et al.
#' @references
#' Coleman, J., Katz, E., & Menzel, H. (1966). Medical innovation: A diffusion
#' study (2nd ed.). New York: Bobbs-Merrill
#'
#' Valente, T. W. (1995). Network models of the diffusion of innovations (2nd ed.).
#' Cresskill N.J.: Hampton Press.
"medInnovations"


#' \code{diffnet} version of the Medical Innovation data
#'
#' A directed dynamic graph with 125 vertices and 18 time periods. The attributes
#' in the graph are static and described in \code{\link{medInnovations}}.
#'
#' @format A \code{\link{diffnet}} class object.
"medInnovationsDiffNet"
