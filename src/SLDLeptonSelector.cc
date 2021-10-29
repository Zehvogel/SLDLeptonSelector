#include "SLDLeptonSelector.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <IMPL/LCCollectionVec.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>
#endif // MARLIN_USE_AIDA

using namespace lcio;
using namespace marlin;

SLDLeptonSelector aSLDLeptonSelector;

SLDLeptonSelector::SLDLeptonSelector() : Processor("SLDLeptonSelector")
{

    // modify processor description
    _description = "SLDLeptonSelector selects leptons from SLD";

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection(LCIO::MCPARTICLE,
                            "MCInColname",
                            "Name of the MCParticle input collection",
                            _mcInColName,
                            std::string("MCParticlesSkimmed"));

    registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "PFOInColName",
                            "Name of the PandoraPFO collection",
                            _pfoInColName,
                            std::string("PandoraPFOs"));

    registerInputCollection(LCIO::LCRELATION,
                            "MCTruthRecoLinkColName",
                            "Name of the MCTruthRecoLink collection",
                            _relInColName,
                            std::string("MCTruthRecoLink"));
    
    registerOutputCollection(LCIO::MCPARTICLE,
                            "MCSLDLeptonsOutColName",
                            "Name of the MCParticle output collection",
                            _mcOutColName,
                            std::string("MCSLDLeptons"));
}

void SLDLeptonSelector::init()
{
    streamlog_out(DEBUG) << "init called" << std::endl;
}

void SLDLeptonSelector::processRunHeader(LCRunHeader *run)
{
    // empty
}

void SLDLeptonSelector::processEvent(LCEvent *evt)
{
    //auto particles = dynamic_cast<LCCollectionVec *>(evt->getCollection(_mcInColName));
    LCCollection *particles = evt->getCollection(_mcInColName);
    if (particles == NULL)
        return;
    LCCollection *MCTruthRecoLink = evt->getCollection(_relInColName);

    auto mcOutCol = new LCCollectionVec(LCIO::MCPARTICLE);
    mcOutCol->setSubset();

    auto nav = std::make_unique<LCRelationNavigator>(MCTruthRecoLink);

    int nMCP = particles->getNumberOfElements();
    // TODO: find a less ugly alternative
    for (int i = 0; i < nMCP; i++)
    {
        MCParticle *p = dynamic_cast<MCParticle *>(particles->getElementAt(i));
        // check if particle decayed
        if (p->getGeneratorStatus() != 2)
            continue;

        // check if particle is of interest i.e. b or c hadron
        int pdg = p->getPDG();
        if (!isBOrCHadron(pdg))
            continue;

        // at least this thing returns a vector for once :)
        // get all stable daughter particles
        auto daughters = p->getDaughters();
        EVENT::MCParticleVec stableLeptons{};
        std::vector<int> stableNeutrinoPDGs{};
        for (const auto& d: daughters) {
            if (d->getGeneratorStatus() != 1)
                continue;
            
            int d_pdg = abs(d->getPDG());
            if (d_pdg == 11 || d_pdg == 13)
                stableLeptons.push_back(d);
            else if (d_pdg == 12 || d_pdg == 14)
                stableNeutrinoPDGs.push_back(d->getPDG());
        } 
        
        for (const auto& l: stableLeptons) {
            for (const int j: stableNeutrinoPDGs) {
                if (abs(j == l->getPDG())) {
                    mcOutCol->addElement(l);
                    break;
                }
            }
        }
    }
    evt->addCollection(mcOutCol, _mcOutColName);
}

bool SLDLeptonSelector::isBOrCHadron(int pdg) {
    bool isCHadron = pdg / 100 == 4 || pdg / 1000 == 4;
    bool isBHadron = pdg / 100 == 5 || pdg / 1000 == 5;
    
    return isBHadron || isCHadron;
}

void SLDLeptonSelector::check(LCEvent *evt)
{
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void SLDLeptonSelector::end()
{

    //   std::cout << "SLDLeptonSelector::end()  " << name()
    // 	    << " processed " << _nEvt << " events in " << _nRun << " runs "
    // 	    << std::endl ;
}
