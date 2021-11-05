#include "SLDLeptonSelector.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/LCRelationNavigator.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/LCRelationImpl.h>

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
    
    registerOutputCollection(LCIO::RECONSTRUCTEDPARTICLE,
                            "PFOOutColName",
                            "Name of the PFO output collection",
                            _pfoOutColName,
                            std::string("SLDLeptonsPFOs"));

    registerInputCollection(LCIO::LCRELATION,
                            "SLDLinkColName",
                            "Name of the SLD link output collection",
                            _relOutColName,
                            std::string("SLDMCRecoLink"));

    registerProcessorParameter("OutputFileName",
                                "Name of the LCIO output file",
                                _outFileName,
                                std::string("test_out.slcio"));
}

void SLDLeptonSelector::init()
{
    streamlog_out(DEBUG) << "init called" << std::endl;
    _lcWriter = LCFactory::getInstance()->createLCWriter();
    _lcWriter->open(_outFileName, LCIO::WRITE_NEW);
}

void SLDLeptonSelector::processRunHeader(LCRunHeader *run)
{
    // empty
}

void SLDLeptonSelector::processEvent(LCEvent *evt)
{
    LCCollection *particles = evt->getCollection(_mcInColName);
    // doesn't work :sad:
    //auto particles = dynamic_cast<MCParticleVec *>(evt->getCollection(_mcInColName));
    if (particles == NULL) {
        streamlog_out(DEBUG) << "particles == NULL" << std::endl;
        return;
    }

    auto mcOutCol = new LCCollectionVec(LCIO::MCPARTICLE);
    mcOutCol->setSubset();

    int nMCP = particles->getNumberOfElements();
    // TODO: find a less ugly alternative
    //for (const auto& p: *particles) {
    for (int i = 0; i < nMCP; i++) {
        MCParticle *p = dynamic_cast<MCParticle *>(particles->getElementAt(i));
        // check if particle decayed
        if (p->getGeneratorStatus() != 2)
            continue;

        // check if particle is of interest i.e. b or c hadron
        int pdg = p->getPDG();
        if (!isBOrCHadron(abs(pdg)))
            continue;

        // get all stable daughter particles
        auto daughters = p->getDaughters();
        MCParticleVec stableLeptons{};
        std::vector<int> stableNeutrinoPDGs{};
        for (const auto &d : daughters)
        {
            if (d->getGeneratorStatus() != 1)
                continue;

            int d_pdg = abs(d->getPDG());
            if (d_pdg == 11 || d_pdg == 13)
                stableLeptons.push_back(d);
            else if (d_pdg == 12 || d_pdg == 14)
                stableNeutrinoPDGs.push_back(d->getPDG());
        }

        for (const auto &l : stableLeptons)
        {
            for (const int j : stableNeutrinoPDGs)
            {
                // streamlog_out(DEBUG) << "[" << j << ", " << l->getPDG() << "]" << std::endl;
                if (abs(j + l->getPDG()) == 1)
                {
                    mcOutCol->addElement(l);
                    break;
                }
            }
        }
    }
    if (mcOutCol->getNumberOfElements() < 1)
        return;
    
    LCCollection *MCTruthRecoLink = evt->getCollection(_relInColName);
    auto nav = std::make_unique<LCRelationNavigator>(MCTruthRecoLink);

    auto recoOutCol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);
    recoOutCol->setSubset();
    
    auto linkOutCol = new LCCollectionVec(LCIO::LCRELATION);

    for (const auto o: *mcOutCol) {
        auto mcp = dynamic_cast<MCParticle *>(o);
        auto recos = nav->getRelatedToObjects(mcp);
        auto weights = nav->getRelatedToWeights(mcp);
        size_t highest = 0;
        for (size_t i = 0; i < recos.size(); i++) {
            auto tweight = int(weights[i]) % 10000;
            if (tweight > 5000) {
                highest = i;
                break;
            } else if (tweight > int(weights[highest]) % 10000) {
                highest = i;
            }
        }
        ReconstructedParticle *reco = NULL;
        if (recos.size() > 0) {
            reco = dynamic_cast<ReconstructedParticle *>(recos[highest]);
            recoOutCol->addElement(reco);
            streamlog_out(DEBUG) << "[mc, rec]: [" << mcp->getPDG() << ", " << reco->getType() << "]" << std::endl;
        }
        // streamlog_out(DEBUG) << "reco: " << reco << std::endl;
        auto rel = new LCRelationImpl(mcp, reco);
        linkOutCol->addElement(rel);
    }

    evt->addCollection(mcOutCol, _mcOutColName);
    evt->addCollection(recoOutCol, _pfoOutColName);
    evt->addCollection(linkOutCol, _relOutColName);
    // TODO: write event to output
    _lcWriter->writeEvent(evt);
}

bool SLDLeptonSelector::isBOrCHadron(int pdg)
{
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
    _lcWriter->close();
    delete _lcWriter;
}
