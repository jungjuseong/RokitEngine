//Copyright (c) 2019 Ultimaker B.V.
//OrganRegenEngine is released under the terms of the AGPLv3 or higher.

#include <assert.h>
#include <cmath>
#include <iomanip>
#include <stdarg.h>

#include "Application.h" //To send layer view data.
#include "ExtruderTrain.h"
#include "gcodeExport.h"
#include "PrintFeature.h"
#include "RetractionConfig.h"
#include "Slice.h"
#include "communication/Communication.h" //To send layer view data.
#include "settings/types/LayerIndex.h"
#include "utils/Date.h"
#include "utils/logoutput.h"
#include "utils/string.h" // MMtoStream, PrecisionedDouble
#include "WipeScriptConfig.h"

namespace cura {

std::string transliterate(const std::string& text)
{
    // For now, just replace all non-ascii characters with '?'.
    // This function can be expaned if we need more complex transliteration.
    std::ostringstream stream;
    for (const char& c : text)
    {
        stream << static_cast<char>((c >= 0) ? c : '?');
    }
    return stream.str();
}

GCodeExport::GCodeExport()
: output_stream(&std::cout)
, currentPosition(0,0,MM2INT(20))
, layer_nr(0)
, relative_extrusion(false)
{
    *output_stream << std::fixed;

    current_e_value = 0;
    current_extruder = 0;
    current_fan_speed = -1;

    total_print_times = std::vector<Duration>(static_cast<unsigned char>(PrintFeatureType::NumPrintFeatureTypes), 0.0);

    currentSpeed = 1;
    current_print_acceleration = -1;
    current_travel_acceleration = -1;
    current_jerk = -1;

    is_z_hopped = 0;
    setFlavor(EGCodeFlavor::INVIVO4D6);
    initial_bed_temp = 0;
    build_volume_temperature = 0;
    machine_heated_build_volume = false;

    fan_number = 0;
    use_extruder_offset_to_offset_coords = false;
    machine_name = "";
    machine_buildplate_type = "";
    relative_extrusion = false;
    new_line = "\n";

    total_bounding_box = AABB3D();

    is_traveling = 1;
    has_first_extruder_setting = false;
}

GCodeExport::~GCodeExport()
{
}

void GCodeExport::preSetup(const size_t start_extruder)
{
    current_extruder = start_extruder;

    const Scene& scene = Application::getInstance().current_slice->scene;
    std::vector<MeshGroup>::iterator mesh_group = scene.current_mesh_group;
    setFlavor(mesh_group->settings.get<EGCodeFlavor>("machine_gcode_flavor"));
    use_extruder_offset_to_offset_coords = mesh_group->settings.get<bool>("machine_use_extruder_offset_to_offset_coords");
    const size_t extruder_count = Application::getInstance().current_slice->scene.extruders.size();

    for (size_t extruder_nr = 0; extruder_nr < extruder_count; extruder_nr++)
    {
        const ExtruderTrain& train = scene.extruders[extruder_nr];
        setFilamentDiameter(extruder_nr, train.settings.get<coord_t>("material_diameter"));

        extruder_attr[extruder_nr].last_retraction_prime_speed = train.settings.get<Velocity>("retraction_prime_speed"); // the alternative would be switch_extruder_prime_speed, but dual extrusion might not even be configured...
        extruder_attr[extruder_nr].fan_number = train.settings.get<size_t>("machine_extruder_cooling_fan_number");
    }

    machine_name = mesh_group->settings.get<std::string>("machine_name");
    machine_buildplate_type = mesh_group->settings.get<std::string>("machine_buildplate_type");

    relative_extrusion = mesh_group->settings.get<bool>("relative_extrusion");

    if (flavor == EGCodeFlavor::BFB)
    {
        new_line = "\r\n";
    }
    else 
    {
        new_line = "\n";
    }

    estimateCalculator.setFirmwareDefaults(mesh_group->settings);
}

void GCodeExport::setInitialAndBuildVolumeTemps(const unsigned int start_extruder_nr)
{
    const Scene& scene = Application::getInstance().current_slice->scene;
    const size_t extruder_count = Application::getInstance().current_slice->scene.extruders.size();
    for (size_t extruder_nr = 0; extruder_nr < extruder_count; extruder_nr++)
    {
        const ExtruderTrain& train = scene.extruders[extruder_nr];

        const Temperature print_temp_0 = train.settings.get<Temperature>("material_print_temperature_layer_0");
        const Temperature print_temp_here = (print_temp_0 != 0)? print_temp_0 : train.settings.get<Temperature>("material_print_temperature");
        const Temperature temp = (extruder_nr == start_extruder_nr)? print_temp_here : train.settings.get<Temperature>("material_standby_temperature");
        setInitialTemp(extruder_nr, temp);
    }

    initial_bed_temp = scene.current_mesh_group->settings.get<Temperature>("material_bed_temperature_layer_0");
    machine_heated_build_volume = scene.current_mesh_group->settings.get<bool>("machine_heated_build_volume");
    build_volume_temperature = machine_heated_build_volume ? scene.current_mesh_group->settings.get<Temperature>("build_volume_temperature") : Temperature(0);
}

void GCodeExport::setInitialTemp(int extruder_nr, double temp)
{
    extruder_attr[extruder_nr].initial_temp = temp;
    if (flavor == EGCodeFlavor::GRIFFIN || flavor == EGCodeFlavor::ULTIGCODE)
    {
        extruder_attr[extruder_nr].currentTemperature = temp;
    }
}

const std::string GCodeExport::flavorToString(const EGCodeFlavor& flavor) const
{
    switch (flavor)
    {
        case EGCodeFlavor::BFB:
            return "BFB";
        case EGCodeFlavor::MACH3:
            return "Mach3";
        case EGCodeFlavor::MAKERBOT:
            return "Makerbot";
        case EGCodeFlavor::ULTIGCODE:
            return "UltiGCode";
        case EGCodeFlavor::GRIFFIN:
            return "Griffin";
        case EGCodeFlavor::REPETIER:
            return "Repetier";
        case EGCodeFlavor::REPRAP:
            return "RepRap";
        case EGCodeFlavor::INVIVO4D6:
            return "Invivo4d6";
        case EGCodeFlavor::MARLIN:
        default:
            return "Marlin";
    }
}

std::string GCodeExport::getFileHeader(const std::vector<bool>& extruder_is_used, const Duration* print_time, const std::vector<double>& filament_used, const std::vector<std::string>& mat_ids)
{
    std::ostringstream prefix;

    const size_t extruder_count = Application::getInstance().current_slice->scene.extruders.size();
    switch (flavor)
    {
    case EGCodeFlavor::GRIFFIN:
        prefix << ";START_OF_HEADER" << new_line;
        prefix << ";HEADER_VERSION:0.1" << new_line;
        prefix << ";FLAVOR:" << flavorToString(flavor) << new_line;
        prefix << ";GENERATOR.NAME:Cura_SteamEngine" << new_line;
        prefix << ";GENERATOR.VERSION:" << VERSION << new_line;
        prefix << ";GENERATOR.BUILD_DATE:" << Date::getDate().toStringDashed() << new_line;
        prefix << ";TARGET_MACHINE.NAME:" << transliterate(machine_name) << new_line;

        for (size_t extr_nr = 0; extr_nr < extruder_count; extr_nr++)
        {
            if (!extruder_is_used[extr_nr])
            {
                continue;
            }
            prefix << ";EXTRUDER_TRAIN." << extr_nr << ".INITIAL_TEMPERATURE:" << extruder_attr[extr_nr].initial_temp << new_line;
            if (filament_used.size() == extruder_count)
            {
                prefix << ";EXTRUDER_TRAIN." << extr_nr << ".MATERIAL.VOLUME_USED:" << static_cast<int>(filament_used[extr_nr]) << new_line;
            }
            if (mat_ids.size() == extruder_count && mat_ids[extr_nr] != "")
            {
                prefix << ";EXTRUDER_TRAIN." << extr_nr << ".MATERIAL.GUID:" << mat_ids[extr_nr] << new_line;
            }
            const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[extr_nr].settings;
            prefix << ";EXTRUDER_TRAIN." << extr_nr << ".NOZZLE.DIAMETER:" << extruder_settings.get<double>("machine_nozzle_size") << new_line;
            prefix << ";EXTRUDER_TRAIN." << extr_nr << ".NOZZLE.NAME:" << extruder_settings.get<std::string>("machine_nozzle_id") << new_line;
        }
        prefix << ";BUILD_PLATE.TYPE:" << machine_buildplate_type << new_line;
        prefix << ";BUILD_PLATE.INITIAL_TEMPERATURE:" << initial_bed_temp << new_line;

        if (machine_heated_build_volume)
        {
            prefix << ";BUILD_VOLUME.TEMPERATURE:" << build_volume_temperature << new_line;
        }

        if (print_time)
        {
            prefix << ";PRINT.TIME:" << static_cast<int>(*print_time) << new_line;
        }

        prefix << ";PRINT.GROUPS:" << Application::getInstance().current_slice->scene.mesh_groups.size() << new_line;

        if (total_bounding_box.min.x > total_bounding_box.max.x) //We haven't encountered any movement (yet). This probably means we're command-line slicing.
        {
            //Put some small default in there.
            total_bounding_box.min = Point3(0, 0, 0);
            total_bounding_box.max = Point3(10, 10, 10);
        }
        prefix << ";PRINT.SIZE.MIN.X:" << INT2MM(total_bounding_box.min.x) << new_line;
        prefix << ";PRINT.SIZE.MIN.Y:" << INT2MM(total_bounding_box.min.y) << new_line;
        prefix << ";PRINT.SIZE.MIN.Z:" << INT2MM(total_bounding_box.min.z) << new_line;
        prefix << ";PRINT.SIZE.MAX.X:" << INT2MM(total_bounding_box.max.x) << new_line;
        prefix << ";PRINT.SIZE.MAX.Y:" << INT2MM(total_bounding_box.max.y) << new_line;
        prefix << ";PRINT.SIZE.MAX.Z:" << INT2MM(total_bounding_box.max.z) << new_line;
        prefix << ";END_OF_HEADER" << new_line;
        break;
    default:
        prefix << ";F/W : 7.7.1.x" << new_line;
        prefix << ";TIME:" << ((print_time)? static_cast<int>(*print_time) : 6666) << new_line;

        prefix << ";Filament used: ";
        if (filament_used.size() > 0)
        {
            for (unsigned i = 0; i < filament_used.size(); ++i)
            {
                if (i > 0)
                {
                    prefix << ", ";
                }
                prefix << filament_used[i] / (1000 * extruder_attr[i].filament_area) << "m";
            }
        }
        else
        {
            prefix << "0m";
        }
        prefix << new_line;
        prefix << ";Layer height: " << Application::getInstance().current_slice->scene.current_mesh_group->settings.get<double>("layer_height") << new_line;
        
        prefix << ";MINX:" << INT2MM(total_bounding_box.min.x) << new_line;
        prefix << ";MINY:" << INT2MM(total_bounding_box.min.y) << new_line;
        prefix << ";MINZ:" << INT2MM(total_bounding_box.min.z) << new_line;
        prefix << ";MAXX:" << INT2MM(total_bounding_box.max.x) << new_line;
        prefix << ";MAXY:" << INT2MM(total_bounding_box.max.y) << new_line;
        prefix << ";MAXZ:" << INT2MM(total_bounding_box.max.z) << new_line;
    }

    return prefix.str();
}


void GCodeExport::setLayerNr(unsigned int layer_nr_) {
    layer_nr = layer_nr_;
}

void GCodeExport::setOutputStream(std::ostream* stream)
{
    output_stream = stream;
    *output_stream << std::fixed;
}

bool GCodeExport::getExtruderIsUsed(const int extruder_nr) const
{
    assert(extruder_nr >= 0);
    assert(extruder_nr < MAX_EXTRUDERS);
    return extruder_attr[extruder_nr].is_used;
}

Point GCodeExport::getGcodePos(const coord_t x, const coord_t y, const int extruder_train) const
{
    if (use_extruder_offset_to_offset_coords)
    {
        const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[extruder_train].settings;
        return Point(x - extruder_settings.get<coord_t>("machine_nozzle_offset_x"), y - extruder_settings.get<coord_t>("machine_nozzle_offset_y"));
    }
    else
    {
        return Point(x, y);
    }
}


void GCodeExport::setFlavor(EGCodeFlavor flavor)
{
    this->flavor = flavor;
    if (flavor == EGCodeFlavor::MACH3)
    {
        for(int n=0; n<MAX_EXTRUDERS; n++)
        {
            extruder_attr[n].extruderCharacter = 'A' + n;
        }
    }
    else
    {
        for(int n=0; n<MAX_EXTRUDERS; n++)
        {
            extruder_attr[n].extruderCharacter = 'E';
        }
    }
    if (flavor == EGCodeFlavor::ULTIGCODE)
    {
        is_volumetric = true;
    }
    else
    {
        is_volumetric = false;
    }
}

EGCodeFlavor GCodeExport::getFlavor() const
{
    return flavor;
}

void GCodeExport::setZ(int z)
{
    current_layer_z = z;
}

void GCodeExport::setFlowRateExtrusionSettings(double max_extrusion_offset, double extrusion_offset_factor)
{
    this->max_extrusion_offset = max_extrusion_offset;
    this->extrusion_offset_factor = extrusion_offset_factor;
}

Point3 GCodeExport::getPosition() const
{
    return currentPosition;
}
Point GCodeExport::getPositionXY() const
{
    return Point(currentPosition.x, currentPosition.y);
}

int GCodeExport::getPositionZ() const
{
    return currentPosition.z;
}

int GCodeExport::getExtruderNr() const
{
    return current_extruder;
}

void GCodeExport::setFilamentDiameter(const size_t extruder, const coord_t diameter)
{
    const double r = INT2MM(diameter) / 2.0;
    const double area = M_PI * r * r;
    extruder_attr[extruder].filament_area = area;
}

double GCodeExport::getCurrentExtrudedVolume() const
{
    double extrusion_amount = current_e_value;
    const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[current_extruder].settings;
    if (!extruder_settings.get<bool>("machine_firmware_retract"))
    { // no E values are changed to perform a retraction
        extrusion_amount -= extruder_attr[current_extruder].retraction_e_amount_at_e_start; // subtract the increment in E which was used for the first unretraction instead of extrusion
        extrusion_amount += extruder_attr[current_extruder].retraction_e_amount_current; // add the decrement in E which the filament is behind on extrusion due to the last retraction
    }
    if (is_volumetric)
    {
        return extrusion_amount;
    }
    else
    {
        return extrusion_amount * extruder_attr[current_extruder].filament_area;
    }
}

double GCodeExport::eToMm(double e)
{
    if (is_volumetric)
    {
        return e / extruder_attr[current_extruder].filament_area;
    }
    else
    {
        return e;
    }
}

double GCodeExport::mm3ToE(double mm3)
{
    if (is_volumetric)
    {
        return mm3;
    }
    else
    {
        return mm3 / extruder_attr[current_extruder].filament_area;
    }
}

double GCodeExport::mmToE(double mm)
{
    if (is_volumetric)
    {
        return mm * extruder_attr[current_extruder].filament_area;
    }
    else
    {
        return mm;
    }
}

double GCodeExport::eToMm3(double e, size_t extruder)
{
    if (is_volumetric)
    {
        return e;
    }
    else
    {
        return e * extruder_attr[extruder].filament_area;
    }
}

double GCodeExport::getTotalFilamentUsed(size_t extruder_nr)
{
    if (extruder_nr == current_extruder)
        return extruder_attr[extruder_nr].totalFilament + getCurrentExtrudedVolume();
    return extruder_attr[extruder_nr].totalFilament;
}

std::vector<Duration> GCodeExport::getTotalPrintTimePerFeature()
{
    return total_print_times;
}

double GCodeExport::getSumTotalPrintTimes()
{
    double sum = 0.0;
    for(double item : getTotalPrintTimePerFeature())
    {
        sum += item;
    }
    return sum;
}

void GCodeExport::resetTotalPrintTimeAndFilament()
{
    for(size_t i = 0; i < total_print_times.size(); i++)
    {
        total_print_times[i] = 0.0;
    }
    for(unsigned int e=0; e<MAX_EXTRUDERS; e++)
    {
        extruder_attr[e].totalFilament = 0.0;
        extruder_attr[e].currentTemperature = 0;
        extruder_attr[e].waited_for_temperature = false;
    }
    current_e_value = 0.0;
    estimateCalculator.reset();
}

void GCodeExport::updateTotalPrintTime()
{
    std::vector<Duration> estimates = estimateCalculator.calculate();
    for(size_t i = 0; i < estimates.size(); i++)
    {
        total_print_times[i] += estimates[i];
    }
    estimateCalculator.reset();
    writeTimeComment(getSumTotalPrintTimes());
}

void GCodeExport::writeExtrudersUsed(const std::vector<bool> extruder_is_used)
{
    *output_stream << ";EXTRUDER_USED:";
    for(size_t i = 0; i < extruder_is_used.size(); i++)
    {
        *output_stream << " " << extruder_is_used[i];
    }
    *output_stream << new_line;
}

void GCodeExport::writeComment(const std::string& unsanitized_comment)
{
    const std::string comment = transliterate(unsanitized_comment);

    *output_stream << ";";
    for (unsigned int i = 0; i < comment.length(); i++)
    {
        if (comment[i] == '\n')
        {
            *output_stream << new_line << ";";
        }
        else
        {
            *output_stream << comment[i];
        }
    }
    *output_stream << new_line;
}

void GCodeExport::writeTimeComment(const Duration time)
{
    *output_stream << ";TIME_ELAPSED:" << time << new_line;
}

void GCodeExport::writeTypeComment(const PrintFeatureType& type)
{
    switch (type)
    {
        case PrintFeatureType::OuterWall:
            *output_stream << ";TYPE:WALL-OUTER" << new_line;
            break;
        case PrintFeatureType::InnerWall:
            *output_stream << ";TYPE:WALL-INNER" << new_line;
            break;
        case PrintFeatureType::Skin:
            *output_stream << ";TYPE:SKIN" << new_line;
            break;
        case PrintFeatureType::Support:
            *output_stream << ";TYPE:SUPPORT" << new_line;
            break;
        case PrintFeatureType::SkirtBrim:
            *output_stream << ";TYPE:SKIRT" << new_line;
            break;
        case PrintFeatureType::Infill:
            *output_stream << ";TYPE:FILL" << new_line;
            break;
        case PrintFeatureType::SupportInfill:
            *output_stream << ";TYPE:SUPPORT" << new_line;
            break;
        case PrintFeatureType::SupportInterface:
            *output_stream << ";TYPE:SUPPORT-INTERFACE" << new_line;
            break;
        case PrintFeatureType::PrimeTower:
            *output_stream << ";TYPE:PRIME-TOWER" << new_line;
            break;
        case PrintFeatureType::MoveCombing:
        case PrintFeatureType::MoveRetraction:
        case PrintFeatureType::NoneType:
        case PrintFeatureType::NumPrintFeatureTypes:
            // do nothing
            break;
    }
}


void GCodeExport::writeLayerComment(const LayerIndex layer_nr)
{
    *output_stream << ";LAYER:" << layer_nr << new_line;
}

void GCodeExport::writeLayerCountComment(const size_t layer_count)
{
    *output_stream << ";LAYER_COUNT:" << layer_count << new_line;
}

void GCodeExport::writeLine(const char* line)
{
    *output_stream << line << new_line;
}

void GCodeExport::writeExtrusionMode(bool set_relative_extrusion_mode)
{
    if (set_relative_extrusion_mode)
    {
        *output_stream << "M83 ;relative extrusion mode" << new_line;
    }
    else 
    {
        if (flavor != EGCodeFlavor::INVIVO4D6)
            *output_stream << "M82 ;absolute extrusion mode" << new_line;
    }
}

void GCodeExport::resetExtrusionValue()
{
    const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[getExtruderNr()].settings;
    const std::string nozzle = extruder_settings.get<std::string>("machine_nozzle_id").c_str();    
    
    if (nozzle.substr(0,3).compare("Ext") == 0 || nozzle.substr(0,3).compare("FFF") == 0) {
        if (!relative_extrusion && has_first_extruder_setting)
        {
            *output_stream << "G92 " << extruder_attr[current_extruder].extruderCharacter << "0" << new_line;
        }
    }
    double current_extruded_volume = getCurrentExtrudedVolume();
    extruder_attr[current_extruder].totalFilament += current_extruded_volume;
    for (double& extruded_volume_at_retraction : extruder_attr[current_extruder].extruded_volume_at_previous_n_retractions)
    { // update the extruded_volume_at_previous_n_retractions only of the current extruder, since other extruders don't extrude the current volume
        extruded_volume_at_retraction -= current_extruded_volume;
    }
    current_e_value = 0.0;
    extruder_attr[current_extruder].retraction_e_amount_at_e_start = extruder_attr[current_extruder].retraction_e_amount_current;
}

void GCodeExport::writeDelay(const Duration& time_amount)
{
    *output_stream << "G4 P" << int(time_amount * 1000) << new_line;
    estimateCalculator.addTime(time_amount);
}

void GCodeExport::writeTravel(const Point& p, const Velocity& speed)
{
    writeTravel(Point3(p.X, p.Y, current_layer_z), speed);
}
void GCodeExport::writeExtrusion(const Point& p, const Velocity& speed, double extrusion_mm3_per_mm, PrintFeatureType feature, bool update_extrusion_offset)
{
    writeExtrusion(Point3(p.X, p.Y, current_layer_z), speed, extrusion_mm3_per_mm, feature, update_extrusion_offset);
}

void GCodeExport::writeTravel(const Point3& p, const Velocity& speed)
{
    if (flavor == EGCodeFlavor::BFB)
    {
        writeMoveBFB(p.x, p.y, p.z + is_z_hopped, speed, 0.0, PrintFeatureType::MoveCombing);
        return;
    }
    writeTravel(p.x, p.y, p.z + is_z_hopped, speed);
}

void GCodeExport::writeExtrusion(const Point3& p, const Velocity& speed, double extrusion_mm3_per_mm, PrintFeatureType feature, bool update_extrusion_offset)
{
    if (flavor == EGCodeFlavor::BFB)
    {
        writeMoveBFB(p.x, p.y, p.z, speed, extrusion_mm3_per_mm, feature);
        return;
    }
    writeExtrusion(p.x, p.y, p.z, speed, extrusion_mm3_per_mm, feature, update_extrusion_offset);
}

void GCodeExport::writeMoveBFB(const int x, const int y, const int z, const Velocity& speed, double extrusion_mm3_per_mm, PrintFeatureType feature)
{
    if (std::isinf(extrusion_mm3_per_mm))
    {
        logError("Extrusion rate is infinite!");
        assert(false && "Infinite extrusion move!");
        std::exit(1);
    }
    if (std::isnan(extrusion_mm3_per_mm))
    {
        logError("Extrusion rate is not a number!");
        assert(false && "NaN extrusion move!");
        std::exit(1);
    }

    double extrusion_per_mm = mm3ToE(extrusion_mm3_per_mm);

    Point gcode_pos = getGcodePos(x,y, current_extruder);

    //For Bits From Bytes machines, we need to handle this completely differently. As they do not use E values but RPM values.
    float fspeed = speed * 60;
    float rpm = extrusion_per_mm * speed * 60;
    const float mm_per_rpm = 4.0; //All BFB machines have 4mm per RPM extrusion.
    rpm /= mm_per_rpm;
    if (rpm > 0)
    {
        if (extruder_attr[current_extruder].retraction_e_amount_current)
        {
            if (currentSpeed != double(rpm))
            {
                //fprintf(f, "; %f e-per-mm %d mm-width %d mm/s\n", extrusion_per_mm, lineWidth, speed);
                //fprintf(f, "M108 S%0.1f\r\n", rpm);
                *output_stream << "M108 S" << PrecisionedDouble{1, rpm} << new_line;
                currentSpeed = double(rpm);
            }
            //Add M101 or M201 to enable the proper extruder.
            *output_stream << "M" << int((current_extruder + 1) * 100 + 1) << new_line;
            extruder_attr[current_extruder].retraction_e_amount_current = 0.0;
        }
        //Fix the speed by the actual RPM we are asking, because of rounding errors we cannot get all RPM values, but we have a lot more resolution in the feedrate value.
        // (Trick copied from KISSlicer, thanks Jonathan)
        fspeed *= (rpm / (roundf(rpm * 100) / 100));

        //Increase the extrusion amount to calculate the amount of filament used.
        Point3 diff = Point3(x,y,z) - getPosition();
        
        current_e_value += extrusion_per_mm * diff.vSizeMM();
    }
    else
    {
        //If we are not extruding, check if we still need to disable the extruder. This causes a retraction due to auto-retraction.
        if (!extruder_attr[current_extruder].retraction_e_amount_current)
        {
            *output_stream << "M103" << new_line;
            extruder_attr[current_extruder].retraction_e_amount_current = 1.0; // 1.0 used as stub; BFB doesn't use the actual retraction amount; it performs retraction on the firmware automatically
        }
    }
    *output_stream << "G1 X" << MMtoStream{gcode_pos.X} << " Y" << MMtoStream{gcode_pos.Y} << " Z" << MMtoStream{z};
    *output_stream << " F" << PrecisionedDouble{1, fspeed} << new_line;
    
    currentPosition = Point3(x, y, z);
    estimateCalculator.plan(TimeEstimateCalculator::Position(INT2MM(currentPosition.x), INT2MM(currentPosition.y), INT2MM(currentPosition.z), eToMm(current_e_value)), speed, feature);
}

void GCodeExport::writeTravel(const coord_t& x, const coord_t& y, const coord_t& z, const Velocity& speed)
{
    if (currentPosition.x == x && currentPosition.y == y && currentPosition.z == z)
    {
        return;
    }

#ifdef ASSERT_INSANE_OUTPUT
    assert(speed < 400 && speed > 1); // normal F values occurring in UM2 gcode (this code should not be compiled for release)
    assert(currentPosition != no_point3);
    assert(Point3(x, y, z) != no_point3);
    assert((Point3(x,y,z) - currentPosition).vSize() < MM2INT(1000)); // no crazy positions (this code should not be compiled for release)
#endif //ASSERT_INSANE_OUTPUT

    const PrintFeatureType travel_move_type = extruder_attr[current_extruder].retraction_e_amount_current ? PrintFeatureType::MoveRetraction : PrintFeatureType::MoveCombing;
    const int display_width = extruder_attr[current_extruder].retraction_e_amount_current ? MM2INT(0.2) : MM2INT(0.1);
    const double layer_height = Application::getInstance().current_slice->scene.current_mesh_group->settings.get<double>("layer_height");
    Application::getInstance().communication->sendLineTo(travel_move_type, Point(x, y), display_width, layer_height, speed);

    const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[getExtruderNr()].settings;
    const std::string nozzle = extruder_settings.get<std::string>("machine_nozzle_id").c_str();    
    
    if (nozzle.substr(0,3).compare("Dis") == 0 || nozzle.substr(0,3).compare("Hot") == 0)
        //if (currentSpeed != speed) 
        //{
        if (is_traveling == 0)
        {
            *output_stream << "M330" << new_line;
            is_traveling = 1;
        }
        //}

    *output_stream << "G0";
    writeFXYZE(speed, x, y, z, current_e_value, travel_move_type);
}

void GCodeExport::writeExtrusion(const int x, const int y, const int z, const Velocity& speed, const double extrusion_mm3_per_mm, const PrintFeatureType& feature, const bool update_extrusion_offset)
{
    if (currentPosition.x == x && currentPosition.y == y && currentPosition.z == z)
    {
        return;
    }

#ifdef ASSERT_INSANE_OUTPUT
    assert(speed < 400 && speed > 1); // normal F values occurring in UM2 gcode (this code should not be compiled for release)
    assert(currentPosition != no_point3);
    assert(Point3(x, y, z) != no_point3);
    assert((Point3(x,y,z) - currentPosition).vSize() < MM2INT(1000)); // no crazy positions (this code should not be compiled for release)
    assert(extrusion_mm3_per_mm >= 0.0);
#endif //ASSERT_INSANE_OUTPUT

    if (std::isinf(extrusion_mm3_per_mm))
    {
        logError("Extrusion rate is infinite!");
        assert(false && "Infinite extrusion move!");
        std::exit(1);
    }

    if (std::isnan(extrusion_mm3_per_mm))
    {
        logError("Extrusion rate is not a number!");
        assert(false && "NaN extrusion move!");
        std::exit(1);
    }

    if (extrusion_mm3_per_mm < 0.0)
    {
        logWarning("Warning! Negative extrusion move!\n");
    }

    double extrusion_per_mm = mm3ToE(extrusion_mm3_per_mm);

    if (is_z_hopped > 0)
    {
        writeZhopEnd();
    }

    Point3 diff = Point3(x,y,z) - currentPosition;

    writeUnretractionAndPrime();

    //flow rate compensation
    double extrusion_offset = 0;
    if (diff.vSizeMM())
    {
        extrusion_offset = speed * extrusion_mm3_per_mm * extrusion_offset_factor;
        if (extrusion_offset > max_extrusion_offset)
        {
            extrusion_offset = max_extrusion_offset;
        }
    }
    // write new value of extrusion_offset, which will be remembered.
    if (update_extrusion_offset && (extrusion_offset != current_e_offset))
    {
        current_e_offset = extrusion_offset;
        *output_stream << ";FLOW_RATE_COMPENSATED_OFFSET = " << current_e_offset << new_line;
    }

    extruder_attr[current_extruder].last_e_value_after_wipe += extrusion_per_mm * diff.vSizeMM();
    double new_e_value = current_e_value + extrusion_per_mm * diff.vSizeMM();
    
    const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[getExtruderNr()].settings;
    const std::string nozzle = extruder_settings.get<std::string>("machine_nozzle_id").c_str();    
    
    if (nozzle.substr(0,3).compare("FFF") != 0 && nozzle.substr(0,3).compare("Ext") != 0)
    {
        //if (currentSpeed != speed || is_traveling == 1)
        if (is_traveling == 1)
        {
            *output_stream << "M301" << new_line;
            is_traveling = 0;
        }
    }
    *output_stream << "G1";
    writeFXYZE(speed, x, y, z, new_e_value, feature);
}

void GCodeExport::writeFXYZE(const Velocity& speed, const int x, const int y, const int z, const double e, const PrintFeatureType& feature)
{
    if (currentSpeed != speed)
    {
        *output_stream << " F" << PrecisionedDouble{1, speed * 60};
        currentSpeed = speed;
    }

    Point gcode_pos = getGcodePos(x, y, current_extruder);
    total_bounding_box.include(Point3(gcode_pos.X, gcode_pos.Y, z));

    *output_stream << " X" << MMtoStream{gcode_pos.X} << " Y" << MMtoStream{gcode_pos.Y};
    if (currentPosition.z != z)
    {
        if (current_extruder == 0)  
            *output_stream << new_line << "G0 Z" << MMtoStream{z};
        else
            *output_stream << new_line << "G0 C" << MMtoStream{z};
    }

    const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[current_extruder].settings;
    std::string nozzle = extruder_settings.get<std::string>("machine_nozzle_id");

    if (e + current_e_offset != current_e_value)
    {
        const double output_e = (relative_extrusion)? e + current_e_offset - current_e_value : e + current_e_offset;
        if (nozzle.substr(0,3).compare("FFF") == 0 || nozzle.substr(0,3).compare("Ext") == 0) 
        {
            *output_stream << " E" << PrecisionedDouble{5, output_e};
        }
    }
    *output_stream << new_line;
    
    currentPosition = Point3(x, y, z);
    current_e_value = e;
    estimateCalculator.plan(TimeEstimateCalculator::Position(INT2MM(x), INT2MM(y), INT2MM(z), eToMm(e)), speed, feature);
}

void GCodeExport::writeUnretractionAndPrime()
{
    const double prime_volume = extruder_attr[current_extruder].prime_volume;
    const double prime_volume_e = mm3ToE(prime_volume);

    const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[current_extruder].settings;
    const std::string nozzle = extruder_settings.get<std::string>("machine_nozzle_id").c_str();    

    const double last_retraction_prime_speed = extruder_attr[current_extruder].last_retraction_prime_speed * 60;

    current_e_value += prime_volume_e;
    if (extruder_attr[current_extruder].retraction_e_amount_current)
    {
        current_e_value += extruder_attr[current_extruder].retraction_e_amount_current;
        const double output_e = (relative_extrusion)? extruder_attr[current_extruder].retraction_e_amount_current + prime_volume_e : current_e_value;
        if (nozzle.substr(0,3).compare("FFF") == 0 || nozzle.substr(0,3).compare("Ext") == 0) 
        {        
            *output_stream << "G1 F" << PrecisionedDouble{1, last_retraction_prime_speed} << " E" << PrecisionedDouble{5, output_e};
        }
        else 
        {
            *output_stream << "M330";
            is_traveling = 1;
        }
         *output_stream << " ;(Back-Retraction)" << new_line;
        
        currentSpeed = extruder_attr[current_extruder].last_retraction_prime_speed;
        estimateCalculator.plan(TimeEstimateCalculator::Position(INT2MM(currentPosition.x), INT2MM(currentPosition.y), INT2MM(currentPosition.z), eToMm(current_e_value)), currentSpeed, PrintFeatureType::MoveRetraction);
    }
    else if (prime_volume != 0.0)
    {
        const double output_e = (relative_extrusion)? prime_volume_e : current_e_value;
        *output_stream << "G1 F" << PrecisionedDouble{1, extruder_attr[current_extruder].last_retraction_prime_speed * 60} << " " << extruder_attr[current_extruder].extruderCharacter;
        *output_stream << PrecisionedDouble{5, output_e} << ";(prime_volume)"  << new_line;
        currentSpeed = extruder_attr[current_extruder].last_retraction_prime_speed;
        estimateCalculator.plan(TimeEstimateCalculator::Position(INT2MM(currentPosition.x), INT2MM(currentPosition.y), INT2MM(currentPosition.z), eToMm(current_e_value)), currentSpeed, PrintFeatureType::NoneType);
    }
    extruder_attr[current_extruder].prime_volume = 0.0;
    
    if (getCurrentExtrudedVolume() > 10000.0) //According to https://github.com/Ultimaker/OrganRegenEngine/issues/14 having more then 21m of extrusion causes inaccuracies. So reset it every 10m, just to be sure.
    {
        resetExtrusionValue();
    }
    if (extruder_attr[current_extruder].retraction_e_amount_current)
    {
        extruder_attr[current_extruder].retraction_e_amount_current = 0.0;
    }
}

void GCodeExport::writeRetraction(const RetractionConfig& config, bool force, bool extruder_switch)
{
    ExtruderTrainAttributes& extr_attr = extruder_attr[current_extruder];

    double old_retraction_e_amount = extr_attr.retraction_e_amount_current;
    double new_retraction_e_amount = mmToE(config.distance);
    double retraction_diff_e_amount = old_retraction_e_amount - new_retraction_e_amount;
    if (std::abs(retraction_diff_e_amount) < 0.000001)
    {
        return;
    }

    { // handle retraction limitation
        double current_extruded_volume = getCurrentExtrudedVolume();
        std::deque<double>& extruded_volume_at_previous_n_retractions = extr_attr.extruded_volume_at_previous_n_retractions;
        while (extruded_volume_at_previous_n_retractions.size() > config.retraction_count_max && !extruded_volume_at_previous_n_retractions.empty()) 
        {
            // extruder switch could have introduced data which falls outside the retraction window
            // also the retraction_count_max could have changed between the last retraction and this
            extruded_volume_at_previous_n_retractions.pop_back();
        }
        if (!force && config.retraction_count_max <= 0)
        {
            return;
        }
        if (!force && extruded_volume_at_previous_n_retractions.size() == config.retraction_count_max
            && current_extruded_volume < extruded_volume_at_previous_n_retractions.back() + config.retraction_extrusion_window * extr_attr.filament_area) 
        {
            return;
        }
        extruded_volume_at_previous_n_retractions.push_front(current_extruded_volume);
        if (extruded_volume_at_previous_n_retractions.size() == config.retraction_count_max + 1) 
        {
            extruded_volume_at_previous_n_retractions.pop_back();
        }
    }

    const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[current_extruder].settings;
    const std::string nozzle = extruder_settings.get<std::string>("machine_nozzle_id").c_str();    

    double speed = ((retraction_diff_e_amount < 0.0) ? config.speed : extr_attr.last_retraction_prime_speed) * 60;
    current_e_value += retraction_diff_e_amount;
    const double output_e = (relative_extrusion)? retraction_diff_e_amount : current_e_value;

    if (has_first_extruder_setting) 
    {
        if (nozzle.substr(0,3).compare("FFF") == 0 || nozzle.substr(0,3).compare("Ext") == 0) {
            *output_stream << "G1 F" << PrecisionedDouble{1, speed} << " E" << PrecisionedDouble{5, output_e};
        }
        else
        {
            *output_stream << "M301";
            is_traveling = 0;
        }
        *output_stream << " ;(Retraction)" << new_line;
    }
    currentSpeed = speed;
    estimateCalculator.plan(TimeEstimateCalculator::Position(INT2MM(currentPosition.x), INT2MM(currentPosition.y), INT2MM(currentPosition.z), eToMm(current_e_value)), currentSpeed, PrintFeatureType::MoveRetraction);
    extr_attr.last_retraction_prime_speed = config.primeSpeed;    

    extr_attr.retraction_e_amount_current = new_retraction_e_amount; // suppose that for UM2 the retraction amount in the firmware is equal to the provided amount
    extr_attr.prime_volume += config.prime_volume;
}

void GCodeExport::writeZhopStart(const coord_t hop_height, Velocity speed/*= 0*/)
{
    if (hop_height > 0)
    {
        if (speed == 0)
        {
            const ExtruderTrain& extruder = Application::getInstance().current_slice->scene.extruders[current_extruder];
            speed = extruder.settings.get<Velocity>("speed_z_hop");
        }
        is_z_hopped = hop_height;
        currentSpeed = speed;
        *output_stream << "G1 F" << PrecisionedDouble{1, speed * 60} << " Z" << MMtoStream{current_layer_z + is_z_hopped} << ";(ENGINE) z hop start" << new_line;
        total_bounding_box.includeZ(current_layer_z + is_z_hopped);
        assert(speed > 0.0 && "Z hop speed should be positive.");
    }
}

void GCodeExport::writeZhopEnd(Velocity speed/*= 0*/)
{
    if (is_z_hopped)
    {
        if (speed == 0)
        {
            const ExtruderTrain& extruder = Application::getInstance().current_slice->scene.extruders[current_extruder];
            speed = extruder.settings.get<Velocity>("speed_z_hop");
        }
        is_z_hopped = 0;
        currentPosition.z = current_layer_z;
        currentSpeed = speed;
        *output_stream << "G1 F" << PrecisionedDouble{1, speed * 60} << " Z" << MMtoStream{current_layer_z} << ";(ENGINE)z hop end" << new_line;
        assert(speed > 0.0 && "Z hop speed should be positive.");
    }
}

void GCodeExport::startExtruder(const size_t new_extruder, const bool from_mesh)
{
    extruder_attr[new_extruder].is_used = true;

    const Settings& currentSettings = Application::getInstance().current_slice->scene.extruders[current_extruder].settings;
    const std::string current_nozzle = currentSettings.get<std::string>("machine_nozzle_id"); 

    const Settings& newSettings = Application::getInstance().current_slice->scene.extruders[new_extruder].settings;
    const std::string new_nozzle = newSettings.get<std::string>("machine_nozzle_id");    

    const Settings& scene_settings = Application::getInstance().current_slice->scene.settings;
    const std::string dish_type = scene_settings.get<std::string>("machine_build_dish_type");

    const std::string LEFT_BED = "G54 G0 X0.0 Y0.0 ;(HOTMELT/EXTRUDER)";
    const std::string RIGHT_BED = "G55 G0 X0.0 Y0.0 ;(ROTARY)";
    const char *A_AXIS_POS[6] = {"0.00", "0.00", "-72.00", "72.00", "144.00", "-144.00"};

    struct StartPoint { 
        std::string x;
        std::string y;
    };
    const std::map<std::string, StartPoint> Jump = {
        {"6", {"-19.50","39.00"}}, 
        {"12", {"-26.00","39.00"}},
        {"24", {"-28.95","48.25"}},
        {"48", {"-32.25","45.15"}},
        {"96", {"-31.50","49.50"}}
    };
    
    const bool is_wellplate = (dish_type.substr(0,10).compare("Well Plate") == 0);
    const bool new_is_extruder = (new_nozzle.substr(0,3).compare("FFF") == 0 || new_nozzle.substr(0,3).compare("Ext") == 0);
    const bool current_is_extruder = (current_nozzle.substr(0,3).compare("FFF") == 0 || current_nozzle.substr(0,3).compare("Ext") == 0);

    if (from_mesh && has_first_extruder_setting == false) 
    {
        // Add HOPPING code when dish type is the WellPlate and tool is the Dispenser
        if (is_wellplate) 
        {
            if (new_extruder == 0)
                return;
                
            std::string well = dish_type.substr(11);
            
            *output_stream << ";HOPPING - "<< "Well No: 0 of " << well << '\n';
            *output_stream << RIGHT_BED << '\n';
            *output_stream << "G0 " << "X" << Jump.at(well).x << " Y" << Jump.at(well).y << '\n';
            *output_stream << "G92 X0.00 Y0.00" << '\n';
            *output_stream << ";END" << '\n';

            *output_stream << ";BODY_START" << '\n';
            *output_stream << ";TOOL_SETUP FROM_MESH: " << new_nozzle.c_str() << " - " << new_extruder << new_line;
            *output_stream << "D" + std::to_string(new_extruder) << '\n';
            *output_stream << "G0 A" << A_AXIS_POS[new_extruder] << " F600" << new_line;
            *output_stream << "G0 B15.0 F300" << new_line;
        } 
        else
        {
            *output_stream << ";BODY_START" << '\n';
            *output_stream << ";TOOL_SETUP FROM_MESH: " << new_nozzle.c_str() << " - " << new_extruder << new_line;
            *output_stream << (new_extruder == 0 ? LEFT_BED : RIGHT_BED) << '\n';
            *output_stream << (new_extruder == 0 ? "D6" : "D" + std::to_string(new_extruder)) << '\n';
            if (new_extruder == 0) {
                //*output_stream << "M301" << new_line;
                //is_traveling = 0;                
            }
            else
            {
                *output_stream << "G0 A" << A_AXIS_POS[new_extruder] << " F600" << new_line;
                *output_stream << "G0 B15.0 F300" << new_line;
            }
        }
        *output_stream << ";END:" << new_line;

        has_first_extruder_setting = true;
        current_extruder = new_extruder;
        
        return;
    }

    if (has_first_extruder_setting && new_extruder != current_extruder) // wouldn't be the case on the very first extruder start if it's extruder 0
    {
        if (is_traveling == 0)
        {
            *output_stream << "M330" << new_line;
            is_traveling = 1;
        }
        *output_stream << ";TOOL_SETUP:" << new_nozzle.c_str() << " - " << new_extruder << new_line;;

        if (new_extruder == 0)
        {
            if (current_is_extruder) {
                if (is_traveling == 0) 
                {
                    *output_stream << "M330" << new_line;
                    is_traveling = 1;
                }
            }
            *output_stream << "D6" << new_line;
            if (new_is_extruder) {
                *output_stream << "M301" << new_line;
                is_traveling = 0;
            }
        }
        else // RIGHT Tools
        {
            *output_stream << "D" << new_extruder << new_line;
            *output_stream << "G0 A" << A_AXIS_POS[new_extruder] << " F600" << new_line;
            *output_stream << "G0 B15.0 F300" << new_line;
        }
        *output_stream << ";END" << new_line;        

        current_extruder = new_extruder;
        has_first_extruder_setting = true;
    }
    assert(getCurrentExtrudedVolume() == 0.0 && "Just after an extruder switch we haven't extruded anything yet!");
    resetExtrusionValue(); // zero the E value on the new extruder, just to be sure

    const std::string start_code = Application::getInstance().current_slice->scene.extruders[new_extruder].settings.get<std::string>("machine_extruder_start_code");

    if(!start_code.empty())
    {
        if (relative_extrusion)
        {
            writeExtrusionMode(false); // ensure absolute extrusion mode is set before the start gcode
        }

        writeCode(start_code.c_str());

        if (relative_extrusion)
        {
            writeExtrusionMode(true); // restore relative extrusion mode
        }
    }

    Application::getInstance().communication->setExtruderForSend(Application::getInstance().current_slice->scene.extruders[new_extruder]);
    Application::getInstance().communication->sendCurrentPosition(getPositionXY());

    //Change the Z position so it gets re-written again. We do not know if the switch code modified the Z position.
    currentPosition.z += 1;

    setExtruderFanNumber(new_extruder);
}

void GCodeExport::switchExtruder(size_t new_extruder, const RetractionConfig& retraction_config_old_extruder, coord_t perform_z_hop /*= 0*/)
{
    if (current_extruder == new_extruder)
    {
        return;
    }

    const Settings& old_extruder_settings = Application::getInstance().current_slice->scene.extruders[current_extruder].settings;
    if(old_extruder_settings.get<bool>("retraction_enable"))
    {
        constexpr bool force = true;
        constexpr bool extruder_switch = true;
        writeRetraction(retraction_config_old_extruder, force, extruder_switch);
    }

    if (perform_z_hop > 0)
    {
        writeZhopStart(perform_z_hop);
    }

    resetExtrusionValue(); // zero the E value on the old extruder, so that the current_e_value is registered on the old extruder

    const std::string end_code = old_extruder_settings.get<std::string>("machine_extruder_end_code");

    if(!end_code.empty())
    {
        if (relative_extrusion)
        {
            writeExtrusionMode(false); // ensure absolute extrusion mode is set before the end gcode
        }

        writeCode(end_code.c_str());

        if (relative_extrusion)
        {
            writeExtrusionMode(true); // restore relative extrusion mode
        }
    }

    startExtruder(new_extruder, false);
}

void GCodeExport::writeCode(const char* str)
{
    *output_stream << str << new_line;
}

void GCodeExport::writePrimeTrain(const Velocity& travel_speed)
{
    if (extruder_attr[current_extruder].is_primed)
    { // extruder is already primed once!
        return;
    }
    const Settings& extruder_settings = Application::getInstance().current_slice->scene.extruders[current_extruder].settings;
    if (extruder_settings.get<bool>("prime_blob_enable"))
    { // only move to prime position if we do a blob/poop
        // ideally the prime position would be respected whether we do a blob or not,
        // but the frontend currently doesn't support a value function of an extruder setting depending on an fdmprinter setting,
        // which is needed to automatically ignore the prime position for the printer when blob is disabled
        Point3 prime_pos(extruder_settings.get<coord_t>("extruder_prime_pos_x"), extruder_settings.get<coord_t>("extruder_prime_pos_y"), extruder_settings.get<coord_t>("extruder_prime_pos_z"));
        if (!extruder_settings.get<bool>("extruder_prime_pos_abs"))
        {
            // currentPosition.z can be already z hopped
            prime_pos += Point3(currentPosition.x, currentPosition.y, current_layer_z);
        }
        writeTravel(prime_pos, travel_speed);
    }

    if (flavor == EGCodeFlavor::GRIFFIN)
    {
        bool should_correct_z = false;
        
        std::string command = "G280";
        if (!extruder_settings.get<bool>("prime_blob_enable"))
        {
            command += " S1";  // use S1 to disable prime blob
            should_correct_z = true;
        }
        *output_stream << command << new_line;

        // There was an issue with the S1 strategy parameter, where it would only change the material-position,
        //   as opposed to 'be a prime-blob maneuvre without actually printing the prime blob', as we assumed here.
        // After a chat, the firmware-team decided to change the S1 strategy behaviour,
        //   but since people don't update their firmware at each opportunity, it was decided to fix it here as well.
        if (should_correct_z)
        {
            // Can't output via 'writeTravel', since if this is needed, the value saved for 'current height' will not be correct.
            // For similar reasons, this isn't written to the front-end via command-socket.
            *output_stream << "G0 Z" << MMtoStream{getPositionZ()} << new_line;
        }
    }
    else
    {
        // there is no prime gcode for other firmware versions...
    }

    extruder_attr[current_extruder].is_primed = true;
}

void GCodeExport::setExtruderFanNumber(int extruder)
{
    if (extruder_attr[extruder].fan_number != fan_number)
    {
        fan_number = extruder_attr[extruder].fan_number;
        current_fan_speed = -1; // ensure fan speed gcode gets output for this fan
    }
}

void GCodeExport::writeFanCommand(double speed)
{
    if (std::abs(current_fan_speed - speed) < 0.1)
    {
        return;
    }
    if (flavor == EGCodeFlavor::INVIVO4D6) {
        return;
    }
    
    if(flavor == EGCodeFlavor::MAKERBOT)
    {
        if(speed >= 50)
        {
            *output_stream << "M126 T0" << new_line; //Makerbot cannot PWM the fan speed...
        }
        else
        {
            *output_stream << "M127 T0" << new_line;
        }
    }
    else if (speed > 0)
    {
        *output_stream << "M106 S" << PrecisionedDouble{1, speed * 255 / 100};
        if (fan_number)
        {
            *output_stream << " P" << fan_number;
        }
        *output_stream << new_line;
    }
    else
    {
        *output_stream << "M107";
        if (fan_number)
        {
            *output_stream << " P" << fan_number;
        }
        *output_stream << new_line;
    }

    current_fan_speed = speed;
}

void GCodeExport::writeTemperatureCommand(const size_t extruder, const Temperature& temperature, const bool wait)
{
    const ExtruderTrain& extruder_train = Application::getInstance().current_slice->scene.extruders[extruder];

    if (flavor == EGCodeFlavor::INVIVO4D6)
        return;

    if (!extruder_train.settings.get<bool>("machine_nozzle_temp_enabled"))
    {
        return;
    }

    if (extruder_train.settings.get<bool>("machine_extruders_share_heater"))
    {
        // extruders share a single heater
        if (extruder != current_extruder)
        {
            // ignore all changes to the non-current extruder
            return;
        }

        // sync all extruders with the change to the current extruder
        const size_t extruder_count = Application::getInstance().current_slice->scene.extruders.size();

        for (size_t extruder_nr = 0; extruder_nr < extruder_count; extruder_nr++)
        {
            if (extruder_nr != extruder)
            {
                extruder_attr[extruder_nr].waited_for_temperature = wait;
                extruder_attr[extruder_nr].currentTemperature = temperature;
            }
        }
    }

    if ((!wait || extruder_attr[extruder].waited_for_temperature) && extruder_attr[extruder].currentTemperature == temperature)
    {
        return;
    }

    if (wait && flavor != EGCodeFlavor::MAKERBOT)
    {
        if (flavor == EGCodeFlavor::INVIVO4D6 || flavor == EGCodeFlavor::MARLIN)
        {
            *output_stream << "M105" << new_line; // get temperatures from the last update, the M109 will not let get the target temperature
        }
        *output_stream << "M109";
        extruder_attr[extruder].waited_for_temperature = true;
    }
    else
    {
        *output_stream << "M104";
        extruder_attr[extruder].waited_for_temperature = false;
    }
    if (extruder != current_extruder)
    {
        *output_stream << " T" << extruder;
    }
#ifdef ASSERT_INSANE_OUTPUT
    assert(temperature >= 0);
#endif // ASSERT_INSANE_OUTPUT
    *output_stream << " S" << PrecisionedDouble{1, temperature} << new_line;
    if (wait && flavor == EGCodeFlavor::MAKERBOT)
    {
        //Makerbot doesn't use M109 for heat-and-wait. Instead, use M104 and then wait using M116.
        *output_stream << "M116" << new_line;
    }
    extruder_attr[extruder].currentTemperature = temperature;
}

void GCodeExport::writeBedTemperatureCommand(const Temperature& temperature, const bool wait)
{
    if (flavor == EGCodeFlavor::INVIVO4D6 || flavor == EGCodeFlavor::ULTIGCODE)
    { // The UM2 family doesn't support temperature commands (they are fixed in the firmware)
        return;
    }

    if (wait)
    {
        if (flavor == EGCodeFlavor::MARLIN)
        {
            *output_stream << "M140 S"; // set the temperature, it will be used as target temperature from M105
            *output_stream << PrecisionedDouble{1, temperature} << new_line;
            *output_stream << "M105" << new_line;
        }
        *output_stream << "M190 S";
    }
    else
        *output_stream << "M140 S";
    *output_stream << PrecisionedDouble{1, temperature} << new_line;
}

void GCodeExport::writeBuildVolumeTemperatureCommand(const Temperature& temperature, const bool wait)
{
    if (flavor == EGCodeFlavor::INVIVO4D6 || flavor == EGCodeFlavor::ULTIGCODE || flavor == EGCodeFlavor::GRIFFIN)
    {
        //Ultimaker printers don't support build volume temperature commands.
        return;
    }
    if (wait)
    {
        *output_stream << "M191 S";
    }
    else
    {
        *output_stream << "M141 S";
    }
    *output_stream << PrecisionedDouble{1, temperature} << new_line;
}

void GCodeExport::writePrintAcceleration(const Acceleration& acceleration)
{
    switch (getFlavor())
    {
        case EGCodeFlavor::REPETIER:
            if (current_print_acceleration != acceleration)
            {
                *output_stream << "M201 X" << PrecisionedDouble{0, acceleration} << " Y" << PrecisionedDouble{0, acceleration} << new_line;
            }
            break;
        case EGCodeFlavor::REPRAP:
            if (current_print_acceleration != acceleration)
            {
                *output_stream << "M204 P" << PrecisionedDouble{0, acceleration} << new_line;
            }
            break;
        default:
            // MARLIN, etc. only have one acceleration for both print and travel
            if (current_print_acceleration != acceleration)
            {
                *output_stream << "M204 S" << PrecisionedDouble{0, acceleration} << new_line;
            }
            break;
    }
    current_print_acceleration = acceleration;
    estimateCalculator.setAcceleration(acceleration);
}

void GCodeExport::writeTravelAcceleration(const Acceleration& acceleration)
{
    switch (getFlavor())
    {
        case EGCodeFlavor::REPETIER:
            if (current_travel_acceleration != acceleration)
            {
                *output_stream << "M202 X" << PrecisionedDouble{0, acceleration} << " Y" << PrecisionedDouble{0, acceleration} << new_line;
            }
            break;
        case EGCodeFlavor::REPRAP:
            if (current_travel_acceleration != acceleration)
            {
                *output_stream << "M204 T" << PrecisionedDouble{0, acceleration} << new_line;
            }
            break;
        default:
            // MARLIN, etc. only have one acceleration for both print and travel
            writePrintAcceleration(acceleration);
            break;
    }
    current_travel_acceleration = acceleration;
    estimateCalculator.setAcceleration(acceleration);
}

void GCodeExport::writeJerk(const Velocity& jerk)
{
    if (flavor == EGCodeFlavor::INVIVO4D6) {
        return;
    }

    if (current_jerk != jerk)
    {
        switch (getFlavor())
        {
            case EGCodeFlavor::REPETIER:
                *output_stream << "M207 X" << PrecisionedDouble{2, jerk} << new_line;
                break;
            case EGCodeFlavor::REPRAP:
                *output_stream << "M566 X" << PrecisionedDouble{2, jerk * 60} << " Y" << PrecisionedDouble{2, jerk * 60} << new_line;
                break;
            default:
                *output_stream << "M205 X" << PrecisionedDouble{2, jerk} << " Y" << PrecisionedDouble{2, jerk} << new_line;
                break;
        }
        current_jerk = jerk;
        estimateCalculator.setMaxXyJerk(jerk);
    }
}

void GCodeExport::finalize(const char* endCode)
{
    writeFanCommand(0);
    writeCode(endCode);
    int64_t print_time = getSumTotalPrintTimes();
    int mat_0 = getTotalFilamentUsed(0);
    log("Print time (s): %d\n", print_time);
    log("Print time (hr|min|s): %dh %dm %ds\n", print_time / 60 / 60, (print_time / 60) % 60, print_time % 60);
    log("Filament (mm^3): %d\n", mat_0);
    for(int n=1; n<MAX_EXTRUDERS; n++)
        if (getTotalFilamentUsed(n) > 0)
            log("Filament%d: %d\n", n + 1, int(getTotalFilamentUsed(n)));
    output_stream->flush();
}

double GCodeExport::getExtrudedVolumeAfterLastWipe(size_t extruder)
{
    return eToMm3(extruder_attr[extruder].last_e_value_after_wipe, extruder);
}

void GCodeExport::ResetLastEValueAfterWipe(size_t extruder)
{
    extruder_attr[extruder].last_e_value_after_wipe = 0;
}

void GCodeExport::insertWipeScript(const WipeScriptConfig& wipe_config)
{
    const Point3 prev_position = currentPosition;
    writeComment("WIPE_SCRIPT_BEGIN");

    if (wipe_config.retraction_enable)
    {
        writeRetraction(wipe_config.retraction_config);
    }

    if (wipe_config.hop_enable)
    {
        writeZhopStart(wipe_config.hop_amount, wipe_config.hop_speed);
    }

    writeTravel(Point(wipe_config.brush_pos_x, currentPosition.y), wipe_config.move_speed);
    for (size_t i = 0; i < wipe_config.repeat_count; ++i)
    {
        coord_t x = currentPosition.x + (i % 2 ? -wipe_config.move_distance : wipe_config.move_distance);
        writeTravel(Point(x, currentPosition.y), wipe_config.move_speed);
    }

    writeTravel(prev_position, wipe_config.move_speed);

    if (wipe_config.hop_enable)
    {
        writeZhopEnd(wipe_config.hop_speed);
    }

    if (wipe_config.retraction_enable)
    {
        writeUnretractionAndPrime();
    }

    if (wipe_config.pause > 0)
    {
        writeDelay(wipe_config.pause);
    }

    writeComment("WIPE_SCRIPT_END");
}

}//namespace cura

